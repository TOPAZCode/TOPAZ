MODULE ModMyRecurrence
use ModProcess
use ModMisc
implicit none


INTERFACE OPERATOR (.Ndot.)
   module procedure FourVecDot
END INTERFACE OPERATOR (.Ndot.)

type :: Insertion
   integer,allocatable :: CacheRef(:)
   integer :: NumGluCur
end type Insertion

integer, parameter :: GluCur_MaxCacheSize=50
integer, parameter :: Dv_Max=8
type :: CachedGluCur
   integer(1) :: rIn
   integer(1) :: rOut
   complex(8) :: cur_g(1:Dv_Max)
   complex(8) :: mom(1:Dv_Max)
end type CachedGluCur

type :: GluonInsertions
   type(Insertion),allocatable :: Ins(:)
   integer :: NumIns
   integer :: CacheSize
   type(CachedGluCur),allocatable:: Cache_GluCur(:)
end type GluonInsertions

type(GluonInsertions),allocatable :: GluInsList_f_2f(:,:)


real(8), parameter :: PropCut = 1.0d-8
integer :: Dv,Ds


integer, parameter :: Cache_PartKey(1:705) = (/ 1,2,3,0,4,0,5,6,0,7,8,0,9,0,10,11,12,0,0,0,13,14,15,0,0,0,0,0,16,0,0,17,18,19,0,20,0,21,22,0,23,24,0,0,25,0,26,27,28,0,0,29,30,0,0,0,0,0,31,32,0,0,33,34,35,0,36,0,0,37,38,0,39,40,0,41,0,42,43,44,0,0,45,0,0,0,0,0,46,47,48,0,0,0,49,50,51,0,52,0,53,54,0,55,56,0,57,0,58,59,60,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,61,62,0,0,63,0,64,0,0,65,66,0,0,0,0,0,0,0,0,0,67,68,0,0,0,0,0,0,69,0,0,70,0,71,0,72,0,0,0,0,0,0,0,0,73,0,74,0,75,0,0,76,0,0,0,0,0,0,77,78,0,0,0,0,0,0,0,0,0,79,80,0,0,81,0,82,0,0,83,84,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,85,86,0,0,87,0,88,0,0,89,90,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,91,92,0,0,0,0,0,0,0,0,0,0,93,0,0,0,94,95,0,0,96,0,0,0,97,0,98,0,0,0,0,0,99,0,0,0,100,0,0,0,0,0,101,0,102,0,0,0,0,103,104,0,0,0,0,0,0,105,0,0,106,0,107,0,108,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,109,110,0,0,0,0,0,0,111,0,0,112,0,113,0,114,0,0,0,0,115,116,0,0,0,0,0,0,0,0,0,0,117,0,0,0,118,119,0,0,120,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,121,0,0,122,123,0,0,0,124,0,0,0,0,0,0,0,0,0,0,125,126,0,0,0,0,127,0,128,0,129,0,0,130,0,0,0,0,0,0,131,132,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,133,0,134,0,135,0,0,136,0,0,0,0,0,0,137,138,0,0,0,0,139,0,140,0,0,0,0,0,141,0,0,0,142,0,0,0,0,0,143,0,144,0,0,0,145,0,0,146,147,0,0,0,148,0,0,0,0,0,0,0,0,0,0,149,150,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,151,152,0,0,153,0,154,0,0,155,156,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,157,158,0,0,159,0,160,0,0,161,162,0,0,0,0,0,0,0,0,0,163,164,0,0,0,0,0,0,165,0,0,166,0,167,0,168,0,0,0,0,0,0,0,0,169,0,170,0,171,0,0,172,0,0,0,0,0,0,173,174,0,0,0,0,0,0,0,0,0,175,176,0,0,177,0,178,0,0,179,180 /)

integer,parameter :: PartComb=180
integer,parameter :: HelComb =16

!DEC$ IF (_Caching_Cur_f_2f==1)
complex(8) :: Cache_cur_f_2f_Ds4(1:4,1:PartComb*HelComb)
complex(8) :: Cache_cur_f_2f_Ds8(1:8,1:PartComb*HelComb)
complex(8) :: Cache_cur_f_2f_Ds16(1:16,1:PartComb*HelComb)        ! memory requirement: 720kByte
logical    :: Cache_cur_f_2f_Ds4_filled(1:PartComb*HelComb)
logical    :: Cache_cur_f_2f_Ds8_filled(1:PartComb*HelComb)
logical    :: Cache_cur_f_2f_Ds16_filled(1:PartComb*HelComb)
private    :: Cache_cur_f_2f_Ds4,Cache_cur_f_2f_Ds8,Cache_cur_f_2f_Ds16,Cache_cur_f_2f_Ds4_filled,Cache_cur_f_2f_Ds8_filled,Cache_cur_f_2f_Ds16_filled
!DEC$ ENDIF

!DEC$ IF (_Caching_Cur_g_2f==1)
complex(8) :: Cache_cur_g_2f_Dv4(1:4,1:PartComb*HelComb)
complex(8) :: Cache_cur_g_2f_Dv6(1:6,1:PartComb*HelComb)
complex(8) :: Cache_cur_g_2f_Dv8(1:8,1:PartComb*HelComb)
logical    :: Cache_cur_g_2f_Dv4_filled(1:PartComb*HelComb)
logical    :: Cache_cur_g_2f_Dv6_filled(1:PartComb*HelComb)
logical    :: Cache_cur_g_2f_Dv8_filled(1:PartComb*HelComb)
private    :: Cache_cur_g_2f_Dv4,Cache_cur_g_2f_Dv6,Cache_cur_g_2f_Dv8,Cache_cur_g_2f_Dv4_filled,Cache_cur_g_2f_Dv6_filled,Cache_cur_g_2f_Dv8_filled
!DEC$ ENDIF


public  :: cur_g,cur_f_2f,cur_g_2f,cur_f_4f,cur_g_4f,cur_f_6f,cur_f_2f_massCT,cur_f_4f_massCT
public  :: cur_f_2fW
public  :: cur_f_2fV
public  :: cur_g_2s,cur_g_ssff,cur_g_ffss,  cur_s_2s,cur_s_ssffss,cur_s_sffsss,cur_s_sssffs,cur_s_4s,cur_s_ssff,cur_s_sffs,  cur_f_ffss,cur_f_fssf,cur_f_fffssf,cur_f_fssfff,  cur_s_2s_massCT
public  :: cur_s_sssffs_CLOSEDLOOPCONTRIB,cur_s_sffsss_CLOSEDLOOPCONTRIB
public  :: SetDim,InitCurrCache
private :: PropCut,Dv,Ds,FourVecDot,Insertion,GluInsList_f_2f,AddGluCur,Cache_PartKey
CONTAINS

! TODO: optimize cycle call for on-shell propagators!
! combine PMom1 and PMom2
! introduce prop cutoff in cur_g(...)



SUBROUTINE CopyParticlePtr(InPointer,OutPointer)
implicit none
type(PtrToParticle), intent(in) :: InPointer
type(PtrToParticle), intent(out):: OutPointer

   OutPointer%PartType => InPointer%PartType
   OutPointer%ExtRef => InPointer%ExtRef
   OutPointer%Mass => InPointer%Mass
   OutPointer%Mass2 => InPointer%Mass2
   OutPointer%Helicity => InPointer%Helicity
   OutPointer%Mom => InPointer%Mom
   OutPointer%Pol => InPointer%Pol
return
END SUBROUTINE






SUBROUTINE SetDim(DvIn,DsIn)
implicit none
integer :: DvIn,DsIn

   Dv = DvIn
   Ds = DsIn
return
END SUBROUTINE


FUNCTION SumMom(Particles,i1,i2)
implicit none
complex(8) :: SumMom(1:Dv)
type(PtrToParticle) :: Particles(:)
integer :: i1,i2,j

   SumMom(1:Dv)= (0d0,0d0)
   do j=i1,i2
      SumMom(1:Dv) = SumMom(1:Dv) + Particles(j)%Mom(1:Dv)
   enddo
END FUNCTION


SUBROUTINE InitCurrCache()
implicit none

!DEC$ IF (_Caching_Cur_f_2f==1)
    Cache_cur_f_2f_Ds4_filled(1:180*16) = .false.
    Cache_cur_f_2f_Ds8_filled(1:180*16) = .false.
    Cache_cur_f_2f_Ds16_filled(1:180*16)= .false.
!DEC$ ENDIF

!DEC$ IF (_Caching_Cur_g_2f==1)
    Cache_cur_g_2f_Dv4_filled(1:180*16) = .false.
    Cache_cur_g_2f_Dv6_filled(1:180*16) = .false.
    Cache_cur_g_2f_Dv8_filled(1:180*16) = .false.

!     Cache_cur_f_2f_hits = 0
!     Cache_cur_f_2f_miss = 0
!DEC$ ENDIF


return
END SUBROUTINE


! FUNCTION cur_g_old(Gluons,NumGlu) result(Res)           ! Gluons(:) does not include the OFF-shell gluon,  however NumGlu is the number of all gluons (ON+OFF-shell)!
! implicit none
! complex(8) :: Res(1:Dv)
! integer :: NumGlu,i
! type(PtrToParticle) :: Gluons(1:)
! complex(8) :: GluMom(1:Dv,1:NumGlu-1)
! complex(8) :: GluPol(1:Dv,1:NumGlu-1)
!
! !DEC$ IF (_DebugCheckMyImpl1==1)
!     if( size(Gluons,dim=1).ne.NumGlu-1 ) print *,"wrong number of gluons in cur_g"
! !DEC$ ENDIF
!    do i=1,NumGlu-1
!     GluMom(1:Dv,i) = Gluons(i)%Mom(1:Dv)
!     GluPol(1:Dv,i) = Gluons(i)%Pol(1:Dv)
!    enddo
!    Res(:) = g(GluPol(1:Dv,:),GluMom(1:Dv,:))
! return
! END FUNCTION






!DEC$ ATTRIBUTES INLINE :: FourVecDot
FUNCTION FourVecDot(p1,p2)
implicit none
complex(8), intent(in) :: p1(1:Dv),p2(1:Dv)
complex(8)  :: FourVecDot
integer :: mu

   FourVecDot = p1(1)*p2(1)
!DEC$ UNROLL
   do mu=2,Dv
      FourVecDot = FourVecDot - p1(mu)*p2(mu)
   enddo
return
END FUNCTION FourVecDot

FUNCTION eval_TripVert(k1,k2,v1,v2)
implicit none
complex(8) :: eval_TripVert(1:Dv)
complex(8) :: k1(1:Dv),k2(1:Dv),v1(1:Dv),v2(1:Dv)
complex(8), parameter :: IOverSqrt2=(0d0,1d0)/dsqrt(2d0)

   eval_TripVert(1:Dv) = IOverSqrt2 * ( (k1(1:Dv)-k2(1:Dv))*(v1.Ndot.v2)  &
                           - 2d0*v1(1:Dv)*(k1.Ndot.v2)  &
                           + 2d0*v2(1:Dv)*(k2.Ndot.v1) )
return
END FUNCTION eval_TripVert

FUNCTION eval_QuadVert(k1,k2,k3)
implicit none
complex(8) :: eval_QuadVert(1:Dv)
complex(8) :: k1(1:Dv),k2(1:Dv),k3(1:Dv)
complex(8), parameter :: I=(0d0,1d0)

   eval_QuadVert(1:Dv) = I * (-k1(1:Dv)*(k2.Ndot.k3)*0.5d0  &
                           + k2(1:Dv)*(k1.Ndot.k3)  &
                           - k3(1:Dv)*(k1.Ndot.k2)*0.5d0 )
return
END FUNCTION eval_QuadVert

!DEC$ ATTRIBUTES INLINE :: linear_map
FUNCTION linear_map(i1,i2,Ngluons)
implicit none
integer :: linear_map,i1,i2,Ngluons

   linear_map = i2+Ngluons*(i2-i1)-((i2-i1)*(i2-i1+1))/2
return
END FUNCTION linear_map



FUNCTION cur_g(Gluons,NumGlu)!  note the off-shell gluon has be counted in NumGlu
implicit none
type(PtrToParticle) :: Gluons(1:)
complex(8) :: cur_g(1:Dv)
integer :: Ngluons,NumGlu
complex(8) :: glu_subcur(1:Dv,1:36)  ! max. 8 gluons allowed
complex(8) :: mom_sum(1:Dv,1:36), PropFactor,PropDenom
integer :: ind0,ind1,ind2,ind3,j,l,mu
integer :: a,b,i1,i2

!DEC$ IF (_DebugWriteCurrents==1)
character :: parts(20)*4
integer :: parti(20)

    do i1=1,NumGlu-1
       if( Gluons(i1)%ExtRef.eq.-1 ) then
 	  exit
       else
          parti(i1)=Gluons(i1)%ExtRef
       endif
       if(i1.eq.NumGlu-1) print parti(1:NumGlu-1)
!       write (parts(1:20), '(I20)') parti(1:NumGlu-1)
    enddo
!DEC$ ENDIF

   Ngluons = NumGlu-1

   do a=0,Ngluons-1
      do b=1,Ngluons-a

         i1 = b
         i2 = a+b
         ind0 = linear_map(i1,i2,Ngluons)

         if (i1.eq.i2) then
            glu_subcur(1:Dv,ind0) = Gluons(i1)%Pol(1:Dv)
            mom_sum(1:Dv,ind0)    = Gluons(i1)%Mom(1:Dv)
         else
            mom_sum(1:Dv,ind0) = mom_sum(1:Dv,ind0-Ngluons+i2-i1-1) + Gluons(i2)%Mom(1:Dv)
            if ( i1 .ne. 1 .or. i2 .ne. Ngluons ) then
               PropDenom = mom_sum(1:Dv,ind0).Ndot.mom_sum(1:Dv,ind0)
               if( abs(PropDenom).lt.PropCut ) cycle
               PropFactor = (0d0,-1d0)/PropDenom
            else
               PropFactor = 1d0
            endif
            do mu=1,Dv
               glu_subcur(mu,ind0) = 0d0
            enddo
            do j=i1,i2-1
               ind1 = linear_map(i1,j,Ngluons)
               ind2 = linear_map(j+1,i2,Ngluons)
               glu_subcur(1:Dv,ind0) = glu_subcur(1:Dv,ind0) +            &
                   eval_TripVert( mom_sum(1:Dv,ind1),mom_sum(1:Dv,ind2),  &
                                   glu_subcur(1:Dv,ind1),glu_subcur(1:Dv,ind2)) * PropFactor
            enddo
            do j=i1,i2-2
               do l=j+1,i2-1
                  ind1 = linear_map(i1,j,Ngluons)
                  ind2 = linear_map(j+1,l,Ngluons)
                  ind3 = linear_map(l+1,i2,Ngluons)
                  glu_subcur(1:Dv,ind0) = glu_subcur(1:Dv,ind0) +  &
                      eval_QuadVert(                                       &
                        glu_subcur(1:Dv,ind1),glu_subcur(1:Dv,ind2),glu_subcur(1:Dv,ind3)) * PropFactor
               enddo
            enddo
         endif
      enddo
   enddo
   cur_g(1:Dv) = glu_subcur(1:Dv,ind0)

return
END FUNCTION




FUNCTION AddGluCur(Cache,CacheCounter,rIn,rOut)
implicit none
integer :: AddGluCur
integer :: i,j
integer(1) :: rIn,rOut
type(CachedGluCur) :: Cache(1:GluCur_MaxCacheSize)
integer :: CacheCounter

   do i=1,CacheCounter
      if( Cache(i)%rIn.ne.rIn .or. Cache(i)%rOut.ne.rOut ) cycle
      AddGluCur=i   ! is only evaluated for a cache hit
      return
   enddo

!  above is only evaluated for a cache miss
   CacheCounter = CacheCounter +1
   if(CacheCounter.gt.GluCur_MaxCacheSize) then
      print *, "Cache list is full"
      AddGluCur=0
   endif
   Cache(CacheCounter)%rIn=rIn
   Cache(CacheCounter)%rOut=rOut
   AddGluCur=CacheCounter
return
END FUNCTION





SUBROUTINE Init_cur_2_2f(NGluMax)
! NGluMax: max. number of gluons on each side
! convention: GluInsList_f_2f(ng1,ng2)%NumIns: number of gluon current insertions
!             GluInsList_f_2f(ng1,ng2)%Cache_GluCur(i): cache line i of type(CachedGluCur)
!             GluInsList_f_2f(ng1,ng2)%Ins(i)%CacheRef(:): insertion i, ordered list of references to cache lines,
!                                                          if>0:left insertion, if<0:right insertion
!             GluInsList_f_2f(ng1,ng2)%Ins(i)%NumGluCur: number of gluon currents for insertion i
!
!               |
!               |-------5
!               |
!        4------|       /6
!               |      /
!               |-----O---7
!               |      \
!               |       \8
!        3------|
!               |
!        2\     |
!          \----|
!        1/     |
!               |
!               X
use ModPermutations
use ModIntegerPartition
#define verbose_ 0
implicit none
integer :: NGluMax,ngLeft,ngRight,NumPart1,NumPart2,NumComb,NumPerm,NumIns,i,j,c1,c2,c3,l,k,n
integer :: LeftCounter,RightCounter,CacheRef,AllocStatus
integer, allocatable :: Partition1(:,:),Partition2(:,:), TheSet(:), PermutSet(:,:)
integer(1), allocatable :: TmpMomList(:)
integer(1) :: rIn,rOut
type(GluonInsertions) :: GluInsertions
type(CachedGluCur) :: Cache_GluCur(1:GluCur_MaxCacheSize)
integer :: CacheCounter

 if( allocated(GluInsList_f_2f) ) print *, "GluInsList_f_2f is already allocated"
 allocate(GluInsList_f_2f(0:NGluMax,0:NGluMax),stat=AllocStatus)
 if( AllocStatus .ne. 0 ) print *, "Memory allocation in TheTree%NumGlu"

 do ngLeft=0,NGluMax
 do ngRight=0,NGluMax

!DEC$ IF( verbose_>=1)
   print *,""
   print *, "Left-Right: ",ngLeft,ngRight
!DEC$ ENDIF

   NumPart1 = PartitionFunction(ngLeft)
   allocate(Partition1(1:NumPart1,1:ngLeft))

   NumPart2 = PartitionFunction(ngRight)
   allocate(Partition2(1:NumPart2,1:ngRight))

!  get Partitions
   call GetIntPartition(ngLeft,Partition1)
   call GetIntPartition(ngRight,Partition2)
   Partition2(1:NumPart2,1:ngRight) = -Partition2(1:NumPart2,1:ngRight)

   NumComb = NumPart1*NumPart2

!  (1) dry run to determine number of insertions
   NumIns=0;
   do i=1,NumPart1
   do j=1,NumPart2

!     dry run to determine number of non-zero entries
      c2=1
      do l=1,ngLeft
         if( Partition1(i,l).ne.0 ) then
	    c2=c2+1
	 endif
      enddo
      do l=1,ngRight
         if( Partition2(j,l).ne.0 ) then
	    c2=c2+1
	 endif
      enddo

!     build insertion set and remove zeros
      allocate(TheSet(1:c2-1))
      c2=1
      do l=1,ngLeft
         if( Partition1(i,l).ne.0 ) then
	    TheSet(c2) = Partition1(i,l)
	    c2=c2+1
	 endif
      enddo
      do l=1,ngRight
         if( Partition2(j,l).ne.0 ) then
	    TheSet(c2) = Partition2(j,l)
	    c2=c2+1
	 endif
      enddo
!     generate permutations of all insertions
      allocate(PermutSet(1:Factorial(c2-1),1:c2-1))
      PermutSet(1:Factorial(c2-1),1:c2-1) = 0
      call GetPermutation(TheSet,PermutSet,.true.,NumPerm)
      NumIns=NumIns+NumPerm

      deallocate(PermutSet)
      deallocate(TheSet)
   enddo
   enddo


!  initialize cache for gluon currents
   CacheCounter=0
   do i=1,GluCur_MaxCacheSize
     Cache_GluCur(i)%rIn =0
     Cache_GluCur(i)%rOut=0
   enddo

!  (2) actual run to build up momentum insertion list
   allocate(GluInsertions%Ins(1:NumIns))
   c1=1;
   do i=1,NumPart1
   do j=1,NumPart2

!     dry run to determine number of non-zero entries
      c2=1
      do l=1,ngLeft
         if( Partition1(i,l).ne.0 ) then
	    c2=c2+1
      endif
      enddo
      do l=1,ngRight
         if( Partition2(j,l).ne.0 ) then
	    c2=c2+1
      endif
      enddo

!     build insertion set and remove zeros
      allocate(TheSet(1:c2-1))
      c2=1
      do l=1,ngLeft
         if( Partition1(i,l).ne.0 ) then
	    TheSet(c2) = Partition1(i,l)
	    c2=c2+1
	 endif
      enddo
      do l=1,ngRight
         if( Partition2(j,l).ne.0 ) then
	    TheSet(c2) = Partition2(j,l)
	    c2=c2+1
	 endif
      enddo

!     generate permutations of all insertions
      allocate(PermutSet(1:Factorial(c2-1),1:c2-1))
      PermutSet(1:Factorial(c2-1),1:c2-1) = 0
      call GetPermutation(TheSet,PermutSet,.true.,NumPerm)
!DEC$ IF( verbose_==2)
      do l=1,NumPerm
         print *, "PermutSet: ",PermutSet(l,1:c2-1)
      enddo
!DEC$ ENDIF

!     generate momentum insertion list for every permutation and save cache reference
      do l=1,NumPerm
      	 allocate(GluInsertions%Ins(c1)%CacheRef(1:c2-1))
         GluInsertions%Ins(c1)%NumGluCur = c2-1
	 allocate(TmpMomList(1:ngLeft+ngRight))

	 c3=1; LeftCounter=1; RightCounter=ngLeft+ngRight;
	 do k=1,c2-1   ! loop over all entries in permut.set, e.g. PermutSet(l,:)=(2,1,-3,1,-1)

      if( PermutSet(l,k).gt.0 ) then
        do n=1,abs(PermutSet(l,k))
		   TmpMomList(c3)=LeftCounter
		   c3=c3+1; LeftCounter=LeftCounter+1;
        enddo
        rIn =TmpMomList(c3-abs(PermutSet(l,k)))
        rOut=TmpMomList(c3-1)
		  GluInsertions%Ins(c1)%CacheRef(k) = +AddGluCur(Cache_GluCur,CacheCounter,rIn,rOut)

      elseif( PermutSet(l,k).lt.0 ) then
        c3=c3+abs(PermutSet(l,k))-1
        do n=1,abs(PermutSet(l,k))
		     TmpMomList(c3)=RightCounter
		     c3=c3-1; RightCounter=RightCounter-1;
        enddo
        c3=c3+abs(PermutSet(l,k))+1
        rIn =TmpMomList(c3-abs(PermutSet(l,k)))
        rOut=TmpMomList(c3-1)
        GluInsertions%Ins(c1)%CacheRef(k) = -AddGluCur(Cache_GluCur,CacheCounter,rIn,rOut)
      endif
	 enddo
	 deallocate(TmpMomList)
	 c1=c1+1
      enddo
      GluInsertions%NumIns = c1-1

      deallocate(PermutSet)
      deallocate(TheSet)
   enddo
   enddo
   deallocate(Partition1)
   deallocate(Partition2)

   allocate(GluInsertions%Cache_GluCur(1:CacheCounter))
   GluInsertions%Cache_GluCur = Cache_GluCur
   GluInsertions%CacheSize = CacheCounter
   GluInsList_f_2f(ngLeft,ngRight) = GluInsertions

   deallocate(GluInsertions%Ins)
   deallocate(GluInsertions%Cache_GluCur)
 enddo
 enddo

return
END SUBROUTINE







FUNCTION cur_f_2f(Gluons,Quarks,Quark1PartType,NumGlu) result(Res)           ! Quarks(:) does not include the OFF-shell quark
implicit none
complex(8) :: Res(1:Ds)
integer :: NumGlu(0:2),i,rIn,rOut,Quark1PartType
type(PtrToParticle) :: Gluons(1:),Quarks(2:2)      ! off-shell quark is not included in Quarks(:)
complex(8) :: GluMom(1:Dv,NumGlu(0)), QuarkMom(1:Dv)
complex(8) :: GluPol(1:Dv,NumGlu(0)), QuarkPol(1:Ds)
character :: FerFla1*3,FerFla2*3
integer :: PartKey,HelKey,CurrKey,Hel_Tmp




!DEC$ IF (_Caching_Cur_f_2f==1)
   if( NumGlu(0).eq.2 .or. NumGlu(0).eq.3 ) then
      PartKey = 0
      HelKey  = 0
      CurrKey =-1
      do i=NumGlu(2),1,-1
         if( Gluons(NumGlu(1)+i)%ExtRef.eq.-1 ) goto 13
         PartKey = PartKey + Gluons(NumGlu(1)+i)%ExtRef   * 5**(NumGlu(2)-i)
         if( Gluons(NumGlu(1)+i)%Helicity.eq.+1 ) then
            Hel_Tmp = 1
         else
            Hel_Tmp = 0
         endif
         HelKey  = HelKey  +  Hel_Tmp * 2**(NumGlu(2)-i)
      enddo
      if( Quarks(2)%ExtRef.eq.-1 ) goto 13
      PartKey = PartKey + Quarks(2)%ExtRef   * 5**(NumGlu(2))
      if( Quarks(2)%Helicity.eq.+1 ) then
         Hel_Tmp = 1
      else
         Hel_Tmp = 0
      endif
      HelKey  = HelKey  +  Hel_Tmp * 2**(NumGlu(2))
      do i=NumGlu(1),1,-1
         if( Gluons(i)%ExtRef.eq.-1 ) goto 13
         PartKey = PartKey + Gluons(i)%ExtRef   * 5**(NumGlu(2)+NumGlu(1)-i+1)
         if( Gluons(i)%Helicity.eq.+1 ) then
            Hel_Tmp = 1
         else
            Hel_Tmp = 0
         endif
         HelKey  = HelKey  + Hel_Tmp * 2**(NumGlu(2)+NumGlu(1)-i+1)
      enddo
      CurrKey = Cache_PartKey(PartKey-37) + 180*HelKey

      if( Ds.eq.4 .and. Cache_cur_f_2f_Ds4_filled(CurrKey) ) then
          Res(1:4) = Cache_cur_f_2f_Ds4(1:4,CurrKey)
!           Cache_cur_f_2f_hits = Cache_cur_f_2f_hits + 1
          return
      elseif( Ds.eq.8 .and. Cache_cur_f_2f_Ds8_filled(CurrKey) ) then
          Res(1:8) = Cache_cur_f_2f_Ds8(1:8,CurrKey)
          return
      elseif( Ds.eq.16 .and. Cache_cur_f_2f_Ds16_filled(CurrKey) ) then
          Res(1:16) = Cache_cur_f_2f_Ds16(1:16,CurrKey)
          return
      endif
!    if( Ds.eq.4) Cache_cur_f_2f_miss = Cache_cur_f_2f_miss + 1

13 continue
   endif
!DEC$ ENDIF



!DEC$ IF (_DebugGeneralChecks==1)
   if( Quarks(2)%PartType .eq.0 .or. .not.IsAQuark(Quarks(2)%PartType)) then
      print *, "Error in cur_f_2f"
      stop
   endif
!DEC$ ENDIF
!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2).ne.0 ) print *, "wrong number of gluons in cur_f_2f"
    if( Quarks(2)%PartType.ne.-Quark1PartType ) print *, "unequal flavors in cur_f_2f"
!DEC$ ENDIF


   do i=1,NumGlu(0)
    GluMom(1:Dv,i) = Gluons(i)%Mom(1:Dv)
    GluPol(1:Dv,i) = Gluons(i)%Pol(1:Dv)
   enddo
   QuarkMom(1:Dv) = Quarks(2)%Mom(1:Dv)
   QuarkPol(1:Ds) = Quarks(2)%Pol(1:Ds)

   if( abs(Quark1PartType).eq.5 ) then
    FerFla1="top"
   elseif( abs(Quark1PartType).eq.6 ) then
    FerFla1="bot"
   elseif( abs(Quark1PartType).eq.3 ) then
    FerFla1="chm"
   else
    FerFla1="str"
   endif

   if( abs(Quarks(2)%PartType).eq.5 ) then
    FerFla2="top"
   elseif( abs(Quarks(2)%PartType).eq.6 ) then
    FerFla2="bot"
   elseif( abs(Quarks(2)%PartType).eq.3 ) then
    FerFla2="chm"
   else
    FerFla2="str"
   endif

   rIn =1
   rOut=NumGlu(0)
   if( Quarks(2)%PartType .gt.0 ) then      !    X----->----
      Res(:) = f(GluPol(1:Dv,rIn:rOut),GluMom(1:Dv,rIn:rOut),QuarkPol(1:Ds),QuarkMom(1:Dv),Quarks(2)%Mass,FerFla2,FerFla1,NumGlu(1))
   else                                     !    X-----<----
      Res(:) = bf(GluPol(1:Dv,rIn:rOut),GluMom(1:Dv,rIn:rOut),QuarkPol(1:Ds),QuarkMom(1:Dv),Quarks(2)%Mass,FerFla2,FerFla1,NumGlu(1))
   endif


!DEC$ IF (_Caching_Cur_f_2f==1)
  if( NumGlu(0).eq.2 .or. NumGlu(0).eq.3 ) then

      if( CurrKey.eq.-1 ) return
      if( Ds.eq.4 ) then
           Cache_cur_f_2f_Ds4(1:4,CurrKey) = Res(1:4)
           Cache_cur_f_2f_Ds4_filled(CurrKey) = .true.
      elseif( Ds.eq.8 ) then
          Cache_cur_f_2f_Ds8(1:8,CurrKey) = Res(1:8)
          Cache_cur_f_2f_Ds8_filled(CurrKey) = .true.
      elseif( Ds.eq.16 ) then
          Cache_cur_f_2f_Ds16(1:16,CurrKey) = Res(1:16)
          Cache_cur_f_2f_Ds16_filled(CurrKey) = .true.
      endif
  endif
!DEC$ ENDIF


return
END FUNCTION





FUNCTION cur_f_2f_new(Gluons,Quarks,Quark1PartType,NumGlu) result(Res)           ! Quarks(:) does not include the OFF-shell quark
implicit none
complex(8) :: Res(1:Ds)
integer :: NumGlu(0:2),CacheLine,Quark1PartType,Top,GluIns,Ref
integer :: rIn,rOut
type(PtrToParticle) :: Gluons(1:),Quarks(2:2)      ! off-shell quark is not included in Quarks(:)
complex(8) :: spi1(1:Ds),PMom1(1:Dv),PropFac1

!DEC$ IF (_DebugGeneralChecks==1)
   if( Quarks(2)%PartType .eq.0 .or. .not.IsAQuark(Quarks(2)%PartType)) then
      print *, "Error in cur_f_2f"
      stop
   endif
!DEC$ ENDIF
!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2).ne.0 ) print *, "wrong number of gluons in cur_f_2f"
    if( Quarks(2)%PartType.ne.-Quark1PartType ) print *, "unequal flavors in cur_f_2f"
!DEC$ ENDIF

   Res(1:Ds)=(0d0,0d0)
   if( NumGlu(0).eq.0 ) then       ! no gluon
      Res(1:Ds) = Quarks(2)%Pol(1:Ds)

   elseif( NumGlu(0).eq.1 .and. NumGlu(1).eq.1) then       ! one gluon
      if( Quarks(2)%PartType.gt.0) then
           Res(1:Ds) = vbqg(Quarks(2)%Pol(1:Ds),Gluons(1)%Pol(1:Dv))
      elseif( Quarks(2)%PartType.lt.0) then
           Res(1:Ds) = vgq(Gluons(1)%Pol(1:Dv),Quarks(2)%Pol(1:Ds))
      endif
   elseif( NumGlu(0).eq.1 .and. NumGlu(2).eq.1) then
      if( Quarks(2)%PartType.gt.0) then
           Res(1:Ds) = vgbq(Gluons(1)%Pol(1:Dv),Quarks(2)%Pol(1:Ds))
      elseif( Quarks(2)%PartType.lt.0) then
           Res(1:Ds) = vqg(Quarks(2)%Pol(1:Ds),Gluons(1)%Pol(1:Dv))
      endif

   else       ! many gluons
!     calc gluon currents in cache
!       print *, "ng1,ng2",NumGlu(1),NumGlu(2)
      do CacheLine=1,GluInsList_f_2f(NumGlu(1),NumGlu(2))%CacheSize
         rIn = GluInsList_f_2f(NumGlu(1),NumGlu(2))%Cache_GluCur(CacheLine)%rIn
         rOut= GluInsList_f_2f(NumGlu(1),NumGlu(2))%Cache_GluCur(CacheLine)%rOut
!          print *, "cache line:",cacheline,":",rin,rout
         if( rIn.ne.rOut ) then
            PMom1(1:Dv)=SumMom(Gluons,rIn,rOut)
            PropFac1=(0d0,-1d0)/sc_(PMom1(1:Dv),PMom1(1:Dv))               ! can be replaced by pointer assignments
            GluInsList_f_2f(NumGlu(1),NumGlu(2))%Cache_GluCur(CacheLine)%Cur_g(1:Dv) = cur_g(Gluons(rIn:rOut),rOut-rIn+2)*PropFac1
            GluInsList_f_2f(NumGlu(1),NumGlu(2))%Cache_GluCur(CacheLine)%Mom(1:Dv)   = PMom1(1:Dv)
         else
            GluInsList_f_2f(NumGlu(1),NumGlu(2))%Cache_GluCur(CacheLine)%Cur_g(1:Dv) = Gluons(rIn)%Pol(1:Dv)
            GluInsList_f_2f(NumGlu(1),NumGlu(2))%Cache_GluCur(CacheLine)%Mom(1:Dv)   = Gluons(rIn)%Mom(1:Dv)
         endif
      enddo

!     calc spinor strings
      do Top=1,GluInsList_f_2f(NumGlu(1),NumGlu(2))%NumIns                       ! loop over topologies
         PMom1(1:Dv) = Quarks(2)%Mom(1:Dv)
         spi1(1:Ds)  = Quarks(2)%Pol(1:Ds)
!          print *, "Top:",Top
         do GluIns=GluInsList_f_2f(NumGlu(1),NumGlu(2))%Ins(Top)%NumGluCur,2,-1
            Ref = GluInsList_f_2f(NumGlu(1),NumGlu(2))%Ins(Top)%CacheRef(GluIns)
!             print *, "Ref1:",Ref
            PMom1(1:Dv) = PMom1(1:Dv) + GluInsList_f_2f(NumGlu(1),NumGlu(2))%Cache_GluCur(abs(Ref))%mom(1:Dv)
            PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(2)%Mass2)
            if( abs(sc_(PMom1,PMom1)-Quarks(2)%Mass2).lt.PropCut ) cycle
            if( Quarks(2)%PartType.gt.0) then
               if( Ref.gt.0 ) then
                  spi1(1:Ds) = vgq(GluInsList_f_2f(NumGlu(1),NumGlu(2))%Cache_GluCur(abs(Ref))%cur_g(1:Dv),spi1(1:Ds))
                  spi1(1:Ds) = (+spb2_(spi1,PMom1)+Quarks(2)%Mass*spi1(1:Ds) )*PropFac1
               else!if( Ref.lt.0 ) then
                  spi1(1:Ds) = vqg(spi1(1:Ds),GluInsList_f_2f(NumGlu(1),NumGlu(2))%Cache_GluCur(abs(Ref))%cur_g(1:Dv))
                  spi1(1:Ds) = (+spb2_(spi1,PMom1)+Quarks(2)%Mass*spi1(1:Ds) )*PropFac1
               endif
            else!if( Quarks(2)%PartType.lt.0) then
               if( Ref.gt.0 ) then
                  spi1(1:Ds) = vgbq(GluInsList_f_2f(NumGlu(1),NumGlu(2))%Cache_GluCur(abs(Ref))%cur_g(1:Dv),spi1(1:Ds))
                  spi1(1:Ds) = (-spi2_(PMom1,spi1)+Quarks(2)%Mass*spi1(1:Ds) )*PropFac1
               else!if( Ref.lt.0 ) then
                  spi1(1:Ds) = vbqg(spi1(1:Ds),GluInsList_f_2f(NumGlu(1),NumGlu(2))%Cache_GluCur(abs(Ref))%cur_g(1:Dv))
                  spi1(1:Ds) = (-spi2_(PMom1,spi1)+Quarks(2)%Mass*spi1(1:Ds) )*PropFac1
               endif
            endif
         enddo
!        last insertion (without propagator)
         Ref = GluInsList_f_2f(NumGlu(1),NumGlu(2))%Ins(Top)%CacheRef(1)
!             print *, "Ref2:",Ref
         if( Quarks(2)%PartType.gt.0) then
            if( Ref.gt.0 ) then
               spi1(1:Ds) = vgq(GluInsList_f_2f(NumGlu(1),NumGlu(2))%Cache_GluCur(abs(Ref))%cur_g(1:Dv),spi1(1:Ds))
            else!if( Ref.lt.0 ) then
               spi1(1:Ds) = vqg(spi1(1:Ds),GluInsList_f_2f(NumGlu(1),NumGlu(2))%Cache_GluCur(abs(Ref))%cur_g(1:Dv))
            endif
         else!if( Quarks(2)%PartType.lt.0) then
            if( Ref.gt.0 ) then
               spi1(1:Ds) = vgbq(GluInsList_f_2f(NumGlu(1),NumGlu(2))%Cache_GluCur(abs(Ref))%cur_g(1:Dv),spi1(1:Ds))
            else!if( Ref.lt.0 ) then
               spi1(1:Ds) = vbqg(spi1(1:Ds),GluInsList_f_2f(NumGlu(1),NumGlu(2))%Cache_GluCur(abs(Ref))%cur_g(1:Dv))
            endif
         endif
         Res(1:Ds) = Res(1:Ds) + spi1(1:Ds)
      enddo

   endif
END FUNCTION



FUNCTION cur_g_2f(Gluons,Quarks,NumGlu) result(Res)           ! Gluons(:) does not include the OFF-shell gluon,  however NumGlu(0) is the number of all gluons
implicit none
integer :: NumGlu(0:3),i,counter
type(PtrToParticle) :: Gluons(1:),Quarks(1:2)
integer :: rIn,rOut,n1a,n1b,n2a,n2b,n3a,n3b
integer,target :: TmpExtRef
complex(8) :: Res(1:Dv)
complex(8) :: u1(1:Ds),ubar2(1:Ds)
complex(8),target :: Eps1(1:Dv)
complex(8) :: Eps2(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(3)+1)
complex(8) :: PMom1(1:Dv),PMom2(1:Dv),PMom4(1:Dv)
complex(8),target :: PMom3(1:Dv)
complex(8) :: PropFac1,PropFac2,PropFac3,PropFac4
integer :: PartKey,HelKey,CurrKey,Hel_Tmp



!DEC$ IF (_Caching_Cur_g_2f==1)
   if( NumGlu(0).eq.3 .or. NumGlu(0).eq.2 ) then
      PartKey = 0
      HelKey  = 0
      CurrKey =-1
      do i=NumGlu(3),1,-1
         if( Gluons(NumGlu(1)+NumGlu(2)+i)%ExtRef.eq.-1 ) goto 14
         PartKey = PartKey + Gluons(NumGlu(1)+NumGlu(2)+i)%ExtRef   * 5**(NumGlu(3)-i)
         if( Gluons(NumGlu(1)+NumGlu(2)+i)%Helicity.eq.+1 ) then
            Hel_Tmp = 1
         else
            Hel_Tmp = 0
         endif
         HelKey  = HelKey  +  Hel_Tmp * 2**(NumGlu(3)-i)
      enddo
      if( Quarks(2)%ExtRef.eq.-1 ) goto 14
      PartKey = PartKey + Quarks(2)%ExtRef   * 5**(NumGlu(3))
      if( Quarks(2)%Helicity.eq.+1 ) then
         Hel_Tmp = 1
      else
         Hel_Tmp = 0
      endif
      HelKey  = HelKey  +  Hel_Tmp * 2**(NumGlu(3))


      do i=NumGlu(2),1,-1
         if( Gluons(NumGlu(1)+i)%ExtRef.eq.-1 ) goto 14
         PartKey = PartKey + Gluons(NumGlu(1)+i)%ExtRef   * 5**(NumGlu(3)+NumGlu(2)-i+1)
         if( Gluons(NumGlu(1)+i)%Helicity.eq.+1 ) then
            Hel_Tmp = 1
         else
            Hel_Tmp = 0
         endif
         HelKey  = HelKey  +  Hel_Tmp * 2**(NumGlu(3)+NumGlu(2)-i+1)
      enddo
      if( Quarks(1)%ExtRef.eq.-1 ) goto 14
      PartKey = PartKey + Quarks(1)%ExtRef   * 5**(NumGlu(3)+NumGlu(2)+1)
      if( Quarks(1)%Helicity.eq.+1 ) then
         Hel_Tmp = 1
      else
         Hel_Tmp = 0
      endif
      HelKey  = HelKey  +  Hel_Tmp * 2**(NumGlu(3)+NumGlu(2)+1)

      do i=NumGlu(1),1,-1
         if( Gluons(i)%ExtRef.eq.-1 ) goto 14
         PartKey = PartKey + Gluons(i)%ExtRef   * 5**(NumGlu(3)+NumGlu(2)+NumGlu(1)-i+2)
         if( Gluons(i)%Helicity.eq.+1 ) then
            Hel_Tmp = 1
         else
            Hel_Tmp = 0
         endif
         HelKey  = HelKey  + Hel_Tmp * 2**(NumGlu(3)+NumGlu(2)+NumGlu(1)-i+2)
      enddo
      CurrKey = Cache_PartKey(PartKey-37) + 180*HelKey


      if( Dv.eq.4 .and. Cache_cur_g_2f_Dv4_filled(CurrKey) ) then
          Res(1:4) = Cache_cur_g_2f_Dv4(1:4,CurrKey)
!           Cache_cur_f_2f_hits = Cache_cur_f_2f_hits + 1
          return
      elseif( Dv.eq.6 .and. Cache_cur_g_2f_Dv6_filled(CurrKey) ) then
          Res(1:6) = Cache_cur_g_2f_Dv6(1:6,CurrKey)
          return
      elseif( Dv.eq.8 .and. Cache_cur_g_2f_Dv8_filled(CurrKey) ) then
          Res(1:8) = Cache_cur_g_2f_Dv8(1:8,CurrKey)
          return
      endif
!    if( Ds.eq.4) Cache_cur_f_2f_miss = Cache_cur_f_2f_miss + 1

14 continue
   endif
!DEC$ ENDIF




!DEC$ IF (_DebugCheckMyImpl1==1)
   if(Quarks(1)%PartType*Quarks(2)%PartType.ge.0) print *,"Error in cur_g_2f: wrong PartTypes"
!DEC$ ENDIF
!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-1-NumGlu(1)-NumGlu(2)-NumGlu(3).ne.0 ) print *, "wrong number of gluons in cur_g_2f"
!DEC$ ENDIF

   res = (0d0,0d0)
   if( Quarks(1)%PartType .ne. -Quarks(2)%PartType ) return
   do n1a=0,NumGlu(1)
   do n3a=0,NumGlu(3)
   do n2a=0,NumGlu(2)
      n1b=NumGlu(1)-n1a
      n2b=NumGlu(2)-n2a
      n3b=NumGlu(3)-n3a
      ! Fer1
      rIn=n1a+1
      rOut=NumGlu(1)+n2a
      PMom1(:) = Quarks(1)%Mom(:)+ SumMom(Gluons,rIn,rOut)
      u1(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(1:1),-Quarks(1)%PartType,(/n1b+n2a,n1b,n2a/))
      if(n1b.ge.1 .or. n2a.ge.1) then
         PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(1)%Mass2)
         if( abs(sc_(PMom1,PMom1)-Quarks(1)%Mass2).lt.PropCut ) cycle
         if( Quarks(1)%PartType.lt.0 ) then
            u1(:) = (-spi2_(PMom1,u1)+Quarks(1)%Mass*u1(:) )*PropFac1
         else
            u1(:) = (+spb2_(u1,PMom1)+Quarks(1)%Mass*u1(:) )*PropFac1
         endif
      endif

      ! Fer2
      rIn=NumGlu(1)+n2a+1
      rOut=NumGlu(1)+NumGlu(2)+n3a
      PMom2(:) = Quarks(2)%Mom(:)+ SumMom(Gluons,rIn,rOut)
      ubar2(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(2:2),-Quarks(2)%PartType,(/n2b+n3a,n2b,n3a/))
      if(n2b.ge.1 .or. n3a.ge.1) then
         PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(2)%Mass2)
         if( abs(sc_(PMom2,PMom2)-Quarks(2)%Mass2).lt.PropCut ) cycle
         if( Quarks(2)%PartType.lt.0 ) then
            ubar2(:) = (-spi2_(PMom2,ubar2)+Quarks(2)%Mass*ubar2(:) )*PropFac2
         else
            ubar2(:) = (+spb2_(ubar2,PMom2)+Quarks(2)%Mass*ubar2(:) )*PropFac2
         endif
      endif

      if( Quarks(1)%PartType.lt.0 ) then
         Eps1(:)= -vbqq(Dv,ubar2,u1)       ! re-checked
      else
         Eps1(:)= +vbqq(Dv,u1,ubar2)       ! re-checked
      endif
      PMom3(:) = Quarks(1)%Mom(:)+Quarks(2)%Mom(:) + SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+n3a)
      counter=1
      rIn =1
      rOut=n1a
      do i=rIn,rOut
         call CopyParticlePtr(Gluons(i),TmpGluons(counter))
         counter=counter+1
      enddo
      TmpGluons(counter)%Mom => PMom3(:)
      TmpGluons(counter)%Pol => Eps1(:)
      TmpExtRef = -1
      TmpGluons(counter)%ExtRef => TmpExtRef
      counter=counter+1
      rIn =NumGlu(1)+NumGlu(2)+n3a+1
      rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)
      do i=rIn,rOut
         call CopyParticlePtr(Gluons(i),TmpGluons(counter))
         counter=counter+1
      enddo
      Eps2(:) = cur_g(TmpGluons(1:counter-1),1+n1a+n3b+1)


      if(n1a.ge.1 .or. n3b.ge.1) then
         PropFac3 = (0d0,-1d0)/sc_(PMom3,PMom3)
         if( abs(sc_(PMom3,PMom3)).lt.PropCut ) cycle
         Eps2(:) = Eps2(:)*PropFac3
      endif

      Res(:) = Res(:) + Eps2(:)
   enddo
   enddo
   enddo


!DEC$ IF (_Caching_Cur_g_2f==1)
   if( NumGlu(0).eq.3 .or. NumGlu(0).eq.2 ) then

      if( CurrKey.eq.-1 ) return
      if( Dv.eq.4 ) then
           Cache_cur_g_2f_Dv4(1:4,CurrKey) = Res(1:4)
           Cache_cur_g_2f_Dv4_filled(CurrKey) = .true.
      elseif( Dv.eq.6 ) then
          Cache_cur_g_2f_Dv6(1:6,CurrKey) = Res(1:6)
          Cache_cur_g_2f_Dv6_filled(CurrKey) = .true.
      elseif( Dv.eq.8 ) then
          Cache_cur_g_2f_Dv8(1:8,CurrKey) = Res(1:8)
          Cache_cur_g_2f_Dv8_filled(CurrKey) = .true.
      endif
  endif
!DEC$ ENDIF


return
END FUNCTION





FUNCTION cur_f_4f(Gluons,Quarks,Quark1PartType,NumGlu,tag_f) result(res)           ! Quarks(:) does not include the OFF-shell quark
implicit none
integer :: NumGlu(0:4),Quark1PartType
type(PtrToParticle) :: Gluons(1:),Quarks(2:4)
integer :: tag_f
integer,target :: TmpExtRef
complex(8) :: res(1:Ds),tmp(1:Ds)
complex(8) :: ubar1(1:Ds)
complex(8),target :: ubar0(1:Ds)
complex(8) :: eps1(1:Dv)
complex(8) :: eps2(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(4)),TmpQuark(1:1)
complex(8) :: PropFac1,PropFac2
complex(8),target :: pmom1(1:Dv)
complex(8) :: pmom2(1:Dv)
integer :: n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b
integer :: rIn,rOut,i,counter


!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4).ne.0 ) print *, "wrong number of gluons in cur_f_4f"
    if(Quarks(3)%PartType.eq.-Quarks(4)%PartType .and. Quark1PartType.ne.-Quarks(2)%PartType ) print *,"wrong flavor in cur_f_4f (1)"
    if(Quarks(2)%PartType.eq.-Quarks(3)%PartType .and. Quark1PartType.ne.-Quarks(4)%PartType ) print *,"wrong flavor in cur_f_4f (2)"
!DEC$ ENDIF


   Res(:)=(0d0,0d0)

   if( Quark1PartType.eq.-Quarks(2)%PartType .and. Quarks(3)%PartType.eq.-Quarks(4)%PartType ) then
!     (I)
      do n2a=0,NumGlu(2)
      do n4a=0,NumGlu(4)
         n2b = NumGlu(2)-n2a
         n4b = NumGlu(4)-n4a

         rIn =NumGlu(1)+n2a+1
         rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
         Eps2 = cur_g_2f(Gluons(rIn:rOut),Quarks(3:4),(/1+n2b+NumGlu(3)+n4a,n2b,NumGlu(3),n4a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom + Quarks(4)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
         Eps2 = Eps2*PropFac1
         do n1a=0,NumGlu(1)
            n1b = NumGlu(1)-n1a
            ! Fer2
            rIn =n1a+1
            rOut=NumGlu(1)+n2a
            ubar1(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(2:2),-Quarks(2)%PartType,(/n2a+n1b,n1b,n2a/) )
            if(n1b.ge.1 .or. n2a.ge.1) then
               PMom2(:) = Quarks(2)%Mom + SumMom(Gluons,rIn,rOut)  ! can be moved outside the n1a-loop
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(2)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Quarks(2)%Mass2).lt.PropCut ) cycle
               if( Quarks(2)%PartType.lt.0 ) then
                  ubar1(:) = (-spi2_(PMom2,ubar1)+Quarks(2)%Mass*ubar1(:))*PropFac2
               else
                  ubar1(:) = (+spb2_(ubar1,PMom2)+Quarks(2)%Mass*ubar1(:))*PropFac2
               endif
            endif
            if( Quarks(2)%PartType.lt.0 ) then
               ubar0(:) = vbqg(ubar1,eps2)       ! re-checked
            else
               ubar0(:) = vqg(ubar1,eps2)        ! re-checked
            endif

            PMom1 = Quarks(2)%Mom+Quarks(3)%Mom+Quarks(4)%Mom+SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a)
            if(n1a.ge.1 .or. n4b.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(2)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Quarks(2)%Mass2).lt.PropCut ) cycle
               if( Quarks(2)%PartType.lt.0 ) then
                  ubar0(:) = (-spi2_(PMom1,ubar0)+Quarks(2)%Mass*ubar0(:))*PropFac1
               else
                  ubar0(:) = (+spb2_(ubar0,PMom1)+Quarks(2)%Mass*ubar0(:))*PropFac1
               endif
            endif

            TmpQuark(1)%Mom  => PMom1(:)
            TmpQuark(1)%Pol  => ubar0(:)
            TmpQuark(1)%Mass => Quarks(2)%Mass
            TmpQuark(1)%Mass2=> Quarks(2)%Mass2
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpQuark(1)%PartType => Quarks(2)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),-TmpQuark(1)%PartType,(/counter-1,n1a,n4b/) )
            Res(:) = Res(:) + tmp(:)
         enddo
      enddo
      enddo
   endif


!    if( (Quark1PartType.eq.-Quarks(4)%PartType .and. Quarks(2)%PartType.eq.-Quarks(3)%PartType) .AND.  &
!        ((abs(Quark1PartType).ne.abs(Quarks(2)%PartType).and.(tag_f.ne.3)) .or. &
!         (abs(Quark1PartType).eq.abs(Quarks(2)%PartType).and.(tag_f.ne.1.and.tag_f.ne.3))) ) then
!    if( (Quark1PartType.eq.-Quarks(4)%PartType .and. Quarks(2)%PartType.eq.-Quarks(3)%PartType) .AND.  &
!        .not.(Quarks(4)%ExtRef.eq.-1 .and. tag_f.eq.1 .and. abs(Quark1PartType).eq.(Quarks(2)%PartType))  ) then
   if( (Quark1PartType.eq.-Quarks(4)%PartType .and. Quarks(2)%PartType.eq.-Quarks(3)%PartType) .AND.  &
       (Quarks(4)%ExtRef.ne.-1 .or. tag_f.ne.1 .or. abs(Quark1PartType).ne.abs(Quarks(2)%PartType)) &
     ) then
!     (II)
      do n1a=0,NumGlu(1)
      do n3b=0,NumGlu(3)
         n1b = NumGlu(1)-n1a
         n3a = NumGlu(3)-n3b
         rIn =n1a+1
         rOut=NumGlu(1)+NumGlu(2)+n3a
         Eps2 = cur_g_2f(Gluons(rIn:rOut),Quarks(2:3),(/1+n1b+NumGlu(2)+n3a,n1b,NumGlu(2),n3a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
         Eps2 = Eps2*PropFac1

         do n4a=0,NumGlu(4)
            n4b = NumGlu(4)-n4a
            ! Fer4
            rIn =NumGlu(1)+NumGlu(2)+n3a+1
            rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
            ubar1(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(4:4),-Quarks(4)%PartType,(/NumGlu(3)+n4a-n3a,n3b,n4a/) )
            if(n3b.ge.1 .or. n4a.ge.1) then
               PMom2(:) = Quarks(4)%Mom + SumMom(Gluons,rIn,rOut)
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(4)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Quarks(4)%Mass2).lt.PropCut ) cycle
               if( Quarks(4)%PartType.lt.0 ) then
                  ubar1(:) = (-spi2_(PMom2,ubar1)+Quarks(4)%Mass*ubar1(:))*PropFac2
               else
                  ubar1(:) = (+spb2_(ubar1,PMom2)+Quarks(4)%Mass*ubar1(:))*PropFac2
               endif
            endif
            if( Quarks(4)%PartType.lt.0 ) then
               ubar0(:) = vgbq(eps2,ubar1)   !! changed from vqg(ubar1,eps2)       ! re-checked
            else
               ubar0(:) = vgq(eps2,ubar1)   !! changed from vbqg(ubar1,eps2)       ! re-checked
            endif

            PMom1 = Quarks(2)%Mom+Quarks(3)%Mom+Quarks(4)%Mom+SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a)
            if(n1a.ge.1 .or. n4b.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(4)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Quarks(4)%Mass2).lt.PropCut ) cycle
               if( Quarks(4)%PartType.lt.0 ) then
                  ubar0(:) = (-spi2_(PMom1,ubar0)+Quarks(4)%Mass*ubar0(:))*PropFac1
               else
                  ubar0(:) = (+spb2_(ubar0,PMom1)+Quarks(4)%Mass*ubar0(:))*PropFac1
               endif
            endif
            TmpQuark(1)%Mom  => PMom1(:)
            TmpQuark(1)%Pol  => ubar0(:)
            TmpQuark(1)%Mass => Quarks(4)%Mass
            TmpQuark(1)%Mass2=> Quarks(4)%Mass2
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpQuark(1)%PartType => Quarks(4)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),-TmpQuark(1)%PartType,(/counter-1,n1a,n4b/) )
            Res(:) = Res(:) + tmp(:)
         enddo
      enddo
      enddo
   endif
return
END FUNCTION







FUNCTION cur_g_4f(Gluons,Quarks,NumGlu) result(res)           ! Gluons(:) does not include the OFF-shell gluon,  however NumGlu is the number of all gluons
implicit none
integer,intent(in) :: NumGlu(0:5)
type(PtrToParticle) :: Gluons(1:),Quarks(1:)
integer :: na,nb,nc,nd,ne,nf,ng,nh,ni,nj,nk
integer :: rIn,rOut
integer :: tag_f,counter,i
complex(8) :: res(Dv)
type(PtrToParticle) :: TmpGluons(1:2+NumGlu(1)+NumGlu(3)+NumGlu(5))
complex(8),target :: TmpMom1(1:Dv),TmpMom2(1:Dv)
integer,target :: TmpExtRef1,TmpExtRef2
complex(8),target :: Eps1(1:Dv)
complex(8),target :: Eps2(1:Dv)
complex(8) :: Eps3(1:Dv)
complex(8) :: u1(1:Ds)
complex(8) :: ubar2(1:Ds)
complex(8) :: PropFac1,PropFac2,PropFac3,PropFac4
complex(8) :: PMom1(1:Dv)
complex(8) :: PMom2(1:Dv)
complex(8) :: PMom3(1:Dv)
complex(8) :: PMom4(1:Dv)

!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-1-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4)-NumGlu(5).ne.0 ) print *, "wrong number of gluons in cur_g_4f"
!DEC$ ENDIF

      res = (0d0,0d0)
      if( (Quarks(1)%PartType.eq.-Quarks(2)%PartType .and. Quarks(3)%PartType.eq.-Quarks(4)%PartType) &
     .OR. (Quarks(1)%PartType.eq.-Quarks(4)%PartType .and. Quarks(2)%PartType.eq.-Quarks(3)%PartType) ) then
!        type (1)
         do na=0,NumGlu(1)
         do nc=0,NumGlu(2)
         do ne=0,NumGlu(5)
            nb=NumGlu(1)-na
            nd=NumGlu(2)-nc
            nf=NumGlu(5)-ne
            rIn = NumGlu(1)+nc+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne
            u1 = cur_f_4f(Gluons(rIn:rOut),Quarks(2:4),Quarks(1)%PartType,(/nd+NumGlu(3)+NumGlu(4)+ne,nd,NumGlu(3),NumGlu(4),ne/),0)
            PMom2  = SumMom(Gluons,rIn,rOut)  + Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom
            PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(1)%Mass2)
            if( abs(sc_(PMom2,PMom2) - Quarks(1)%Mass2).lt.PropCut ) cycle
            if( Quarks(1)%PartType.lt.0 ) then
               u1 = ( spb2_(u1,PMom2) + Quarks(1)%Mass*u1 )*PropFac2
            else
               u1 = (-spi2_(PMom2,u1) + Quarks(1)%Mass*u1 )*PropFac2
            endif

            rIn = na+1
            rOut= NumGlu(1)+nc
            ubar2 = cur_f_2f(Gluons(rIn:rOut),Quarks(1:1),-Quarks(1)%PartType,(/nb+nc,nb,nc/))
            PMom3 = SumMom(Gluons,rIn,rOut) + Quarks(1)%Mom
            if( nb.ge.1 .or. nc.ge.1 ) then
               PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Quarks(1)%Mass2)
               if( abs(sc_(PMom3,PMom3) - Quarks(1)%Mass2).lt.PropCut ) cycle
               if( Quarks(1)%PartType.lt.0 ) then
                  ubar2 = (-spi2_(PMom3,ubar2) + Quarks(1)%Mass*ubar2 )*PropFac3
               else
                  ubar2 = (+spb2_(ubar2,PMom3) + Quarks(1)%Mass*ubar2 )*PropFac3
               endif
            endif

            if( Quarks(1)%PartType.lt.0 ) then
                Eps1 = -vbqq(Dv,u1,ubar2)       ! re-checked
            else
                Eps1 = +vbqq(Dv,ubar2,u1)       ! re-checked
            endif

            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpMom1(:) = PMom2(:)+PMom3(:)
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps1(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps2(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+1)

            if( na.ge.1 .or. nf.ge.1 ) then
               PropFac1 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
               if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
               Eps2 = Eps2*PropFac1
            endif

            Res = Res + Eps2
         enddo
         enddo
         enddo


!        type (2)
         do na=0,NumGlu(1)
         do nc=0,NumGlu(4)
         do ne=0,NumGlu(5)        ! can be replaced by above ne-loop
            nb=NumGlu(1)-na
            nd=NumGlu(4)-nc
            nf=NumGlu(5)-ne

            rIn = na+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+nc
            u1 = cur_f_4f(Gluons(rIn:rOut),Quarks(1:3),Quarks(4)%PartType,(/nb+NumGlu(2)+NumGlu(3)+nc,nb,NumGlu(2),NumGlu(3),nc/),0)
            PMom2 =  SumMom(Gluons,rIn,rOut) + Quarks(1)%Mom + Quarks(2)%Mom + Quarks(3)%Mom
            PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
            if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
            if( Quarks(4)%PartType.lt.0 ) then
              u1 = (+spb2_(u1,PMom2) + Quarks(4)%Mass*u1 )*PropFac2
            else
              u1 = (-spi2_(PMom2,u1) + Quarks(4)%Mass*u1 )*PropFac2
            endif
            rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+nc+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne
            ubar2 = cur_f_2f(Gluons(rIn:rOut),Quarks(4:4),-Quarks(4)%PartType,(/nd+ne,nd,ne/))
            PMom3  = SumMom(Gluons,rIn,rOut) + Quarks(4)%Mom
            if( nd.ge.1 .or. ne.ge.1 ) then
               PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Quarks(4)%Mass2)
               if( abs(sc_(PMom3,PMom3) - Quarks(4)%Mass2).lt.PropCut ) cycle
               if( Quarks(4)%PartType.lt.0 ) then
                  ubar2 = (-spi2_(PMom3,ubar2) + Quarks(4)%Mass*ubar2 )*PropFac3
               else
                  ubar2 = (+spb2_(ubar2,PMom3) + Quarks(4)%Mass*ubar2 )*PropFac3
               endif
            endif

            if( Quarks(4)%PartType.lt.0 ) then
              Eps1 = +vbqq(Dv,u1,ubar2)       ! re-checked
            else
              Eps1 = -vbqq(Dv,ubar2,u1)       ! re-checked
            endif

            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpMom1(:) = PMom2(:)+PMom3(:)
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps1(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps2(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+1)

            if( na.ge.1 .or. nf.ge.1 ) then
               PropFac1 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
               if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
               Eps2 = Eps2*PropFac1
            endif

            Res = Res + Eps2
         enddo
         enddo
         enddo
      endif



      if( Quarks(1)%PartType.eq.-Quarks(2)%PartType .and. Quarks(3)%PartType.eq.-Quarks(4)%PartType) then
!        type(3)
         do na=0,NumGlu(1)
         do nc=0,NumGlu(2)
         do ne=0,NumGlu(3)
         do nf=0,NumGlu(3)-ne
         do nh=0,NumGlu(4)   ! this loop can be placed after Eps1 has been calculated
         do nj=0,NumGlu(5)
            nb=NumGlu(1)-na
            nd=NumGlu(2)-nc
            ng=NumGlu(3)-ne-nf
            ni=NumGlu(4)-nh
            nk=NumGlu(5)-nj

            rIn = na+1
            rOut= NumGlu(1)+nc
            ubar2 = cur_f_2f(Gluons(rIn:rOut),Quarks(1:1),-Quarks(1)%PartType,(/nb+nc,nb,nc/))
            PMom1 = SumMom(Gluons,rIn,rOut) + Quarks(1)%Mom
            if( nb.ge.1 .or. nc.ge.1 ) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1) - Quarks(1)%Mass2)
               if( abs(sc_(PMom1,PMom1) - Quarks(1)%Mass2).lt.PropCut ) cycle
               if( Quarks(1)%PartType.lt.0 ) then
                  ubar2 = (-spi2_(PMom1,ubar2) + Quarks(1)%Mass*ubar2 )*PropFac1
               else
                  ubar2 = (+spb2_(ubar2,PMom1) + Quarks(1)%Mass*ubar2 )*PropFac1
               endif
            endif
            rIn = NumGlu(1)+nc+1
            rOut= NumGlu(1)+NumGlu(2)+ne
            u1 = cur_f_2f(Gluons(rIn:rOut),Quarks(2:2),-Quarks(2)%PartType,(/nd+ne,nd,ne/))
            PMom2 = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom
            if( nd.ge.1 .or. ne.ge.1 ) then
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
               if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
               if( Quarks(2)%PartType.lt.0 ) then
                  u1 = (-spi2_(PMom2,u1) + Quarks(2)%Mass*u1 )*PropFac2
               else
                  u1 = (+spb2_(u1,PMom2) + Quarks(2)%Mass*u1 )*PropFac2
               endif
            endif

            if( Quarks(2)%PartType.lt.0 ) then
              Eps1 = +vbqq(Dv,ubar2,u1)       ! re-checked
            else
              Eps1 = -vbqq(Dv,u1,ubar2)       ! re-checked
            endif
            TmpMom1 = PMom1 + PMom2
            PropFac3 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
            if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
            Eps1 = Eps1*PropFac3


            rIn = NumGlu(1)+NumGlu(2)+ne+nf+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+nh
            u1 = cur_f_2f(Gluons(rIn:rOut),Quarks(3:3),-Quarks(3)%PartType,(/ng+nh,ng,nh/))
            PMom3 = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom
            if( ng.ge.1 .or. nh.ge.1 ) then
               PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Quarks(3)%Mass2)
               if( abs(sc_(PMom3,PMom3) - Quarks(3)%Mass2).lt.PropCut ) cycle
               if( Quarks(3)%PartType.lt.0 ) then
                  u1 = (-spi2_(PMom3,u1) + Quarks(3)%Mass*u1 )*PropFac3
               else
                  u1 = (+spb2_(u1,PMom3) + Quarks(3)%Mass*u1 )*PropFac3
               endif
            endif
            rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+nh+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+nj
            ubar2 = cur_f_2f(Gluons(rIn:rOut),Quarks(4:4),-Quarks(4)%PartType,(/ni+nj,ni,nj/))
            PMom4 = SumMom(Gluons,rIn,rOut) + Quarks(4)%Mom
            if( ni.ge.1 .or. nj.ge.1 ) then
               PropFac4 = (0d0,1d0)/(sc_(PMom4,PMom4) - Quarks(4)%Mass2)
               if( abs(sc_(PMom4,PMom4) - Quarks(4)%Mass2).lt.PropCut ) cycle
               if( Quarks(4)%PartType.lt.0 ) then
                  ubar2 = (-spi2_(PMom4,ubar2) + Quarks(4)%Mass*ubar2 )*PropFac4
               else
                  ubar2 = (+spb2_(ubar2,PMom4) + Quarks(4)%Mass*ubar2 )*PropFac4
               endif
            endif

            if( Quarks(4)%PartType.lt.0 ) then
               Eps2 = +vbqq(Dv,u1,ubar2)       ! re-checked
            else
               Eps2 = -vbqq(Dv,ubar2,u1)       ! re-checked
            endif
            TmpMom2 = PMom3 + PMom4
            PropFac1 = (0d0,-1d0)/sc_(TmpMom2,TmpMom2)
            if( abs(sc_(TmpMom2,TmpMom2)).lt.PropCut ) cycle
            Eps2 = Eps2*PropFac1


            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps1(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+ne+1
            rOut=NumGlu(1)+NumGlu(2)+ne+nf
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpGluons(counter)%Mom => TmpMom2(:)
            TmpGluons(counter)%Pol => Eps2(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+nj+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps3(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+nk+2)
            Res = Res + Eps3
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
      endif


return
END FUNCTION




FUNCTION cur_f_6f(Gluons,Quarks,Quark1PartType,NumGlu,tag_f) result(res)           ! Quarks(:) does not include the OFF-shell quark
implicit none
integer :: NumGlu(0:6),Quark1PartType,tag_f
type(PtrToParticle) :: Gluons(1:),Quarks(2:6)
integer,target :: TmpPartType,TmpExtRef
complex(8) :: Res(1:Ds),tmp(1:Ds)
! complex(8) :: Res1(1:Ds),Res2(1:Ds),Res3(1:Ds),Res4(1:Ds)
complex(8) :: u1(1:Ds),ubar1(1:Ds)
complex(8),target :: ubar0(1:Ds)
complex(8) :: eps1(1:Dv)
complex(8) :: eps2(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(6)),TmpQuark(1:1)
complex(8) :: PropFac1,PropFac2
complex(8),target :: PMom1(1:Dv)
complex(8),target :: PMom2(1:Dv)
integer :: n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,n5a,n5b,n6a,n6b
integer :: rIn,rOut,i,counter


!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4)-NumGlu(5)-NumGlu(6).ne.0 ) print *, "wrong number of gluons in cur_f_6f"
!DEC$ ENDIF

    Res = (0d0,0d0)
!     Res1=(0d0,0d0); Res2=(0d0,0d0); Res3=(0d0,0d0); Res4=(0d0,0d0);


!   (A)
    if( Quark1PartType.eq.-Quarks(2)%PartType .AND. (Quarks(3)%PartType.eq.-Quarks(4)%PartType .or. Quarks(3)%PartType.eq.-Quarks(6)%PartType) &
        .AND. (Quarks(2)%ExtRef.ne.-1 .or. tag_f.ne.1) &
      ) then

      do n2a=0,NumGlu(2)
      do n6a=0,NumGlu(6)
         n2b = NumGlu(2)-n2a
         n6b = NumGlu(6)-n6a

         rIn =NumGlu(1)+n2a+1
         rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
         Eps2 = cur_g_4f(Gluons(rIn:rOut),Quarks(3:6),(/1+n2b+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a,n2b,NumGlu(3),NumGlu(4),NumGlu(5),n6a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom + Quarks(4)%Mom + Quarks(5)%Mom + Quarks(6)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut) cycle
         Eps2 = Eps2*PropFac1
         do n1a=0,NumGlu(1)
            n1b = NumGlu(1)-n1a
            ! Fer2
            rIn =n1a+1
            rOut=NumGlu(1)+n2a
            ubar1(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(2:2),-Quarks(2)%PartType,(/n1b+n2a,n1b,n2a/) )
            if(n1b.ge.1 .or. n2a.ge.1) then
               PMom2(:) = Quarks(2)%Mom + SumMom(Gluons,rIn,rOut)
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(2)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Quarks(2)%Mass2).lt.PropCut ) cycle
               if( Quarks(2)%PartType.lt.0 ) then
                  ubar1(:) = (-spi2_(PMom2,ubar1)+Quarks(2)%Mass*ubar1(:))*PropFac2
               else
                  ubar1(:) = (+spb2_(ubar1,PMom2)+Quarks(2)%Mass*ubar1(:))*PropFac2
               endif
            endif
            if( Quarks(2)%PartType.lt.0 ) then
               ubar0(:) = vbqg(ubar1,eps2)       ! re-checked
            else
               ubar0(:) = vqg(ubar1,eps2)       ! re-checked
            endif

            rIn = n1a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            PMom1 = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom + Quarks(5)%Mom + Quarks(6)%Mom  ! can be simplified with PMom1(:)
            if(n1a.ge.1 .or. n6b.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(2)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Quarks(2)%Mass2).lt.PropCut ) cycle
               if( Quarks(2)%PartType.lt.0 ) then
                  ubar0(:) = (-spi2_(PMom1,ubar0)+Quarks(2)%Mass*ubar0(:))*PropFac1
               else
                  ubar0(:) = (+spb2_(ubar0,PMom1)+Quarks(2)%Mass*ubar0(:))*PropFac1
               endif
            endif

            TmpQuark(1)%Mom  => PMom1(:)
            TmpQuark(1)%Pol  => ubar0(:)
            TmpQuark(1)%Mass => Quarks(2)%Mass
            TmpQuark(1)%Mass2=> Quarks(2)%Mass2
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpQuark(1)%PartType => Quarks(2)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),-TmpQuark(1)%PartType,(/counter-1,n1a,n6b/) )

            Res(:) = Res(:) + tmp(:)
!             Res1(:) = Res1(:) + tmp(:)
!             print *, "1",tmp(:)
         enddo
      enddo
      enddo
    endif



!   (B)
    if( Quark1PartType.eq.-Quarks(6)%PartType .AND. (Quarks(2)%PartType.eq.-Quarks(5)%PartType .or. Quarks(2)%PartType.eq.-Quarks(3)%PartType) &
        .AND. (Quarks(6)%ExtRef.ne.-1 .or. tag_f.ne.1) &
      ) then
      do n1a=0,NumGlu(1)
      do n5b=0,NumGlu(5)
         n1b = NumGlu(1)-n1a
         n5a = NumGlu(5)-n5b

         rIn =n1a+1
         rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a
         Eps2 = cur_g_4f(Gluons(rIn:rOut),Quarks(2:5),(/1+n1b+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a,n1b,NumGlu(2),NumGlu(3),NumGlu(4),n5a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom + Quarks(5)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
         Eps2 = Eps2*PropFac1

         do n6a=0,NumGlu(6)
            n6b = NumGlu(6)-n6a
            ! Fer6
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a+1
            rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            ubar1(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(6:6),-Quarks(6)%PartType,(/n5b+n6a,n5b,n6a/) )
            if(n5b.ge.1 .or. n6a.ge.1) then
               PMom2(:) = SumMom(Gluons,rIn,rOut) + Quarks(6)%Mom
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(6)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Quarks(6)%Mass2).lt.PropCut ) cycle
               if( Quarks(6)%PartType.lt.0 ) then
                  ubar1(:) = (-spi2_(PMom2,ubar1)+Quarks(6)%Mass*ubar1(:))*PropFac2
               else
                  ubar1(:) = (+spb2_(ubar1,PMom2)+Quarks(6)%Mass*ubar1(:))*PropFac2
               endif
            endif
            if( Quarks(6)%PartType.lt.0 ) then
               ubar0(:) = vgbq(Eps2,ubar1)       ! re-checked
            else
               ubar0(:) = vgq(Eps2,ubar1)       ! re-checked
            endif

            rIn = n1a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            PMom1 = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom + Quarks(5)%Mom + Quarks(6)%Mom
            if(n1a.ge.1 .or. n6b.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(6)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Quarks(6)%Mass2).lt.PropCut ) cycle
               if( Quarks(6)%PartType.lt.0 ) then
                  ubar0(:) = (-spi2_(PMom1,ubar0)+Quarks(6)%Mass*ubar0(:))*PropFac1
               else
                  ubar0(:) = (+spb2_(ubar0,PMom1)+Quarks(6)%Mass*ubar0(:))*PropFac1
               endif
            endif
            TmpQuark(1)%Mom => PMom1(:)
            TmpQuark(1)%Pol => ubar0(:)
            TmpQuark(1)%Mass => Quarks(6)%Mass
            TmpQuark(1)%Mass2=> Quarks(6)%Mass2
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpQuark(1)%PartType => Quarks(6)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),Quark1PartType,(/counter-1,n1a,n6b/) )

            Res(:) = Res(:) + tmp(:)
!             Res2(:) = Res2(:) + tmp(:)
!             print *, "2",tmp(:)
         enddo
      enddo
      enddo
    endif




!   (C)
!     if( Quarks(5)%PartType.eq.-Quarks(6)%PartType .AND. (Quark1PartType.eq.-Quarks(2)%PartType .or. Quark1PartType.eq.-Quarks(4)%PartType) &
!         .AND. .not.(Quarks(6)%ExtRef.eq.-1 .and. tag_f.eq.1) &
!       ) then
!     if( Quark1PartType.eq.-Quarks(2)%PartType .and. Quarks(5)%PartType.eq.-Quarks(6)%PartType .and..not.(Quarks(2)%ExtRef.eq.-1 .and. tag_f.eq.1) &
!    .OR. Quark1PartType.eq.-Quarks(4)%PartType .and. Quarks(5)%PartType.eq.-Quarks(6)%PartType .and..not.(Quarks(4)%ExtRef.eq.-1 .and. tag_f.eq.1) &
!       ) then

    if( Quarks(5)%PartType.eq.-Quarks(6)%PartType .AND. &
        ((Quark1PartType.eq.-Quarks(2)%PartType .and. (Quarks(2)%ExtRef.ne.-1.or.tag_f.ne.1) ) &
    .OR. (Quark1PartType.eq.-Quarks(4)%PartType .and. (Quarks(4)%ExtRef.ne.-1.or.tag_f.ne.1) ))&
      ) then

      do n1a=0,NumGlu(1)
      do n4a=0,NumGlu(4)
      do n6a=0,NumGlu(6)
         n1b = NumGlu(1)-n1a
         n4b = NumGlu(4)-n4a
         n6b = NumGlu(6)-n6a

            rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            Eps2 = cur_g_2f(Gluons(rIn:rOut),Quarks(5:6),(/1+n4b+NumGlu(5)+n6a,n4b,NumGlu(5),n6a/))
            PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom
            PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
            if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
            Eps2 = Eps2*PropFac1

            rIn = n1a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
            u1 = cur_f_4f(Gluons(rIn:rOut),Quarks(2:4),Quark1PartType,(/n1b+NumGlu(2)+NumGlu(3)+n4a,n1b,NumGlu(2),NumGlu(3),n4a/),0)
            PMom2  = SumMom(Gluons,rIn,rOut)  + Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom

            if( Quark1PartType.eq.-Quarks(2)%PartType) then
                if( Quarks(2)%PartType.lt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
                  u1 = (-spi2_(PMom2,u1) + Quarks(2)%Mass*u1 )*PropFac2
                  ubar0 = vbqg(u1,eps2)       ! re-checked
                  rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                  rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom   ! this PMom2 will be re-used below
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
                      ubar0 = (-spi2_(PMom2,ubar0) + Quarks(2)%Mass*ubar0 )*PropFac2
                  endif
                elseif( Quarks(2)%PartType.gt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
                  u1 = ( spb2_(u1,PMom2) + Quarks(2)%Mass*u1 )*PropFac2
                  ubar0 = vqg(u1,eps2)       ! re-checked
                  rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                  rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
                      ubar0 = ( spb2_(ubar0,PMom2) + Quarks(2)%Mass*ubar0 )*PropFac2
                  endif
                endif
                TmpQuark(1)%Mom  => PMom2(:)
                TmpQuark(1)%Pol  => ubar0(:)
                TmpQuark(1)%Mass => Quarks(2)%Mass
                TmpQuark(1)%Mass2=> Quarks(2)%Mass2
            elseif( Quark1PartType.eq.-Quarks(4)%PartType) then
                if( Quarks(4)%PartType.lt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = (-spi2_(PMom2,u1) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vbqg(u1,eps2)       ! re-checked
                  rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                  rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = (-spi2_(PMom2,ubar0) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                elseif( Quarks(4)%PartType.gt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = ( spb2_(u1,PMom2) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vqg(u1,eps2)       ! re-checked
                  rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                  rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = ( spb2_(ubar0,PMom2) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                endif
                TmpQuark(1)%Mom  => PMom2(:)
                TmpQuark(1)%Pol  => ubar0(:)
                TmpQuark(1)%Mass => Quarks(4)%Mass
                TmpQuark(1)%Mass2=> Quarks(4)%Mass2
            endif
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpPartType = -Quark1PartType
            TmpQuark(1)%PartType => TmpPartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),Quark1PartType,(/counter-1,n1a,n6b/) )

            Res(:) = Res(:) + tmp(:)
!             Res3(:) = Res3(:) + tmp(:)
!             print *, "3",tmp(:)
      enddo
      enddo
      enddo
    endif


!   (D)
    if( Quarks(2)%PartType.eq.-Quarks(3)%PartType .AND. ( &
        (Quark1PartType.eq.-Quarks(4)%PartType .and. (Quarks(4)%ExtRef.ne.-1.or.tag_f.ne.1)) &
   .OR. (Quark1PartType.eq.-Quarks(6)%PartType .and. (Quarks(6)%ExtRef.ne.-1.or.tag_f.ne.1))  )) then
!     if( Quark1PartType.eq.-Quarks(4)%PartType .and. Quarks(2)%PartType.eq.-Quarks(3)%PartType .and..not.(Quarks(4)%ExtRef.eq.-1 .and. tag_f.eq.1) &
!     .OR.Quark1PartType.eq.-Quarks(6)%PartType .and. Quarks(2)%PartType.eq.-Quarks(3)%PartType .and..not.(Quarks(6)%ExtRef.eq.-1 .and. tag_f.eq.1) &
!       ) then

      do n1a=0,NumGlu(1)
      do n3a=0,NumGlu(3)
      do n6a=0,NumGlu(6)
         n1b = NumGlu(1)-n1a
         n3b = NumGlu(3)-n3a
         n6b = NumGlu(6)-n6a

            counter=1
            rIn = n1a+1
            rOut= NumGlu(1)+NumGlu(2)+n3a
            Eps2 = cur_g_2f(Gluons(rIn:rOut),Quarks(2:3),(/1+n1b+NumGlu(2)+n3a,n1b,NumGlu(2),n3a/))
            PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
            PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
            if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
            Eps2 = Eps2*PropFac1

            rIn = NumGlu(1)+NumGlu(2)+n3a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            u1 = cur_f_4f(Gluons(rIn:rOut),Quarks(4:6),Quark1PartType,(/n3b+NumGlu(4)+NumGlu(5)+n6a,n3b,NumGlu(4),NumGlu(5),n6a/),tag_f)
            PMom2  = SumMom(Gluons,rIn,rOut)  + Quarks(4)%Mom + Quarks(5)%Mom + Quarks(6)%Mom

            if( Quark1PartType.eq.-Quarks(4)%PartType ) then
                if( Quarks(4)%PartType.lt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = (-spi2_(PMom2,u1) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vgbq(eps2,u1)       ! re-checked
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom   !this PMom2 will be re-used below  ! CAN be written as PMom2=PMom2+PMom1
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = (-spi2_(PMom2,ubar0) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                elseif( Quarks(4)%PartType.gt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = ( spb2_(u1,PMom2) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vgq(eps2,u1)       ! re-checked
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = ( spb2_(ubar0,PMom2) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                endif
                TmpQuark(1)%Mom  => PMom2(:)
                TmpQuark(1)%Pol  => ubar0(:)
                TmpQuark(1)%Mass => Quarks(4)%Mass
                TmpQuark(1)%Mass2=> Quarks(4)%Mass2
            elseif(Quark1PartType.eq.-Quarks(6)%PartType) then
                if( Quarks(6)%PartType.lt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(6)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(6)%Mass2).lt.PropCut ) cycle
                  u1 = (-spi2_(PMom2,u1) + Quarks(6)%Mass*u1 )*PropFac2
                  ubar0 = vgbq(eps2,u1)       ! re-checked
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(6)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(6)%Mass2).lt.PropCut ) cycle
                      ubar0 = (-spi2_(PMom2,ubar0) + Quarks(6)%Mass*ubar0 )*PropFac2
                  endif
                elseif( Quarks(6)%PartType.gt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(6)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(6)%Mass2).lt.PropCut ) cycle
                  u1 = ( spb2_(u1,PMom2) + Quarks(6)%Mass*u1 )*PropFac2
                  ubar0 = vgq(eps2,u1)       ! re-checked
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(6)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(6)%Mass2).lt.PropCut ) cycle
                      ubar0 = ( spb2_(ubar0,PMom2) + Quarks(6)%Mass*ubar0 )*PropFac2
                  endif
                endif
                TmpQuark(1)%Mom  => PMom2(:)
                TmpQuark(1)%Pol  => ubar0(:)
                TmpQuark(1)%Mass => Quarks(6)%Mass
                TmpQuark(1)%Mass2=> Quarks(6)%Mass2
            endif
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpPartType = -Quark1PartType
            TmpQuark(1)%PartType => TmpPartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),Quark1PartType,(/counter-1,n1a,n6b/))

            Res(:)  = Res(:) + Tmp(:)
!             Res4(:) = Res4(:) + Tmp(:)
!             print *, "4",tmp(:)
      enddo
      enddo
      enddo
    endif
!     Res(:) = Res1(:)+Res2(:)+Res3(:)+Res4(:)
!     print *, "Res1:",Res1(1:4)
!     print *, "Res2:",Res2(1:4)
!     print *, "Res3:",Res3(1:4)
!     print *, "Res4:",Res4(1:4)
return
END FUNCTION









FUNCTION cur_f_2f_massCT(Gluons,Quarks,Quark1PartType,NumGlu) result(Res)
implicit none
complex(8) :: Res(1:Ds)
integer :: NumGlu(0:2),i,rIn,rOut,Quark1PartType
type(PtrToParticle) :: Gluons(1:),Quarks(2:2)      ! off-shell quark is not included in Quarks(:)
complex(8) :: GluMom(1:Dv,NumGlu(0)), QuarkMom(1:Dv)
complex(8) :: GluPol(1:Dv,NumGlu(0)), QuarkPol(1:Ds)
character :: FerFla1*3,FerFla2*3


   do i=1,NumGlu(0)
    GluMom(1:Dv,i) = Gluons(i)%Mom(1:Dv)
    GluPol(1:Dv,i) = Gluons(i)%Pol(1:Dv)
   enddo
   QuarkMom(1:Dv) = Quarks(2)%Mom(1:Dv)
   QuarkPol(1:Ds) = Quarks(2)%Pol(1:Ds)

   if( abs(Quark1PartType).eq.5 ) then
    FerFla1="top"
   elseif( abs(Quark1PartType).eq.6 ) then
    FerFla1="bot"
   elseif( abs(Quark1PartType).eq.3 ) then
    FerFla1="chm"
   else
    FerFla1="str"
   endif

   if( abs(Quarks(2)%PartType).eq.5 ) then
    FerFla2="top"
   elseif( abs(Quarks(2)%PartType).eq.6 ) then
    FerFla2="bot"
   elseif( abs(Quarks(2)%PartType).eq.3 ) then
    FerFla2="chm"
   else
    FerFla2="str"
   endif

   rIn =1
   rOut=NumGlu(0)
   if( Quarks(2)%PartType.gt.0 ) then      !    X----->----
      Res(:) = fmCT(GluPol(1:Dv,rIn:rOut),GluMom(1:Dv,rIn:rOut),QuarkPol(1:Ds),QuarkMom(1:Dv),Quarks(2)%Mass,FerFla2,FerFla1,NumGlu(1),.true.)
!       Res(:) = f(GluPol(1:Dv,rIn:rOut),GluMom(1:Dv,rIn:rOut),QuarkPol(1:Ds),QuarkMom(1:Dv),FerFla2,FerFla1,NumGlu(1))
   else                                     !    X-----<----
      print *,"Error in cur_f_2f_CT: Quarks(1)%PartType.gt.0 not implemented"
   endif

return
END FUNCTION







FUNCTION cur_f_4f_massCT(Gluons,Quarks,Quark1PartType,NumGlu) result(res)
implicit none
integer :: NumGlu(0:4),Quark1PartType
type(PtrToParticle) :: Gluons(1:),Quarks(2:4)            ! off-shell quark is not included in Quarks(.)
integer,target :: TmpExtRef,TmpPartType
complex(8) :: res(1:Ds),tmp(1:Ds)
complex(8),target :: spiDr(1:Ds),spiUDr(1:Ds)
complex(8) :: eps1(1:Dv)
complex(8) :: eps2(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(4)),TmpQuark(1:1)
complex(8) :: PropFac1,PropFac2
complex(8),target :: pmom1(1:Dv)
complex(8) :: pmom2(1:Dv)
integer :: n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b
integer :: rIn,rOut,i,counter




   Res(:)=(0d0,0d0)
   if( Quark1PartType.eq.-Quarks(2)%PartType .and. Quarks(3)%PartType.eq.-Quarks(4)%PartType ) then
      if( NumGlu(1).eq.0 .and. NumGlu(2).eq.0 .and. NumGlu(4).eq.0 ) return
!DEC$ IF (_DebugGeneralChecks==1)
      if( Quarks(3)%Mass.ne.0d0 .or. Quarks(4)%Mass.ne.0d0 ) print *, "CT insertions for massive quarks 3,4 not implemented"
!DEC$ ENDIF

      do n2a=0,NumGlu(2)
      do n4a=0,NumGlu(4)
         n2b = NumGlu(2)-n2a
         n4b = NumGlu(4)-n4a

         rIn =NumGlu(1)+n2a+1
         rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
         Eps2 = cur_g_2f(Gluons(rIn:rOut),Quarks(3:4),(/1+n2b+NumGlu(3)+n4a,n2b,NumGlu(3),n4a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom + Quarks(4)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         Eps2 = Eps2*PropFac1

         do n1a=0,NumGlu(1)
            n1b = NumGlu(1)-n1a
            ! Fer2
            rIn =n1a+1
            rOut=NumGlu(1)+n2a
            spiUDr(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(2:2),-Quarks(2)%PartType,(/n2a+n1b,n1b,n2a/))          ! undressed current
            spiDr(:) = (0d0,0d0)
            if(n1b.ge.1 .or. n2a.ge.1) then
               spiDr(:)  = cur_f_2f_massCT(Gluons(rIn:rOut),Quarks(2:2),-Quarks(2)%PartType,(/n2a+n1b,n1b,n2a/))   ! current dressed with mass CT's
               PMom2(:)  = Quarks(2)%Mom + SumMom(Gluons,rIn,rOut)
               PropFac2  = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(2)%Mass2)
               spiDr(:)  = ( spb2_(spiDr,PMom2)+Quarks(2)%Mass*spiDr(:) )*PropFac2
               spiDr(:)  = spiDr(:) + ( spb2_(spiUDr,PMom2)*(2d0*Quarks(2)%Mass) + (sc_(PMom2,PMom2)+Quarks(2)%Mass2)*spiUDr )*PropFac2**2
               spiUDr(:) = ( spb2_(spiUDr,PMom2)+Quarks(2)%Mass*spiUDr(:) )*PropFac2
               spiDr(:)  = vqg(spiDr, eps2)
            endif
            spiUDr(:) = vqg(spiUDr,eps2)

            PMom1 = Quarks(2)%Mom+Quarks(3)%Mom+Quarks(4)%Mom+SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a)
            if(n1a.ge.1 .or. n4b.ge.1) then
               PropFac1  = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(2)%Mass2)
               spiDr(:)  = ( spb2_(spiDr,PMom1)+Quarks(2)%Mass*spiDr(:) )*PropFac1
               spiDr(:)  = spiDr(:) + ( spb2_(spiUDr,PMom1)*(2d0*Quarks(2)%Mass) + (sc_(PMom1,PMom1)+Quarks(2)%Mass2)*spiUDr )*PropFac1**2
               spiUDr(:) = ( spb2_(spiUDr,PMom1)+Quarks(2)%Mass*spiUDr(:) )*PropFac1
            endif

            TmpExtRef = -1
            TmpQuark(1)%Mom  => PMom1(:)
            TmpQuark(1)%Pol  => spiDr(:)
            TmpQuark(1)%Mass => Quarks(2)%Mass
            TmpQuark(1)%Mass2=> Quarks(2)%Mass2
            TmpQuark(1)%ExtRef  => TmpExtRef
            TmpQuark(1)%PartType=> Quarks(2)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo

            tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),-TmpQuark(1)%PartType,(/counter-1,n1a,n4b/) )
            Res(:) = Res(:) + tmp(:)
            if( n1a.ge.1 .or. n4b.ge.1 ) then
              TmpQuark(1)%Pol  => spiUDr(:)
              tmp(:) = cur_f_2f_massCT(TmpGluons(1:counter-1),TmpQuark(1:1),-TmpQuark(1)%PartType,(/counter-1,n1a,n4b/) )
              Res(:) = Res(:) + tmp(:)
            endif
         enddo
      enddo
      enddo




   elseif( Quark1PartType.eq.-Quarks(4)%PartType .and. Quarks(2)%PartType.eq.-Quarks(3)%PartType ) then
      if( NumGlu(1).eq.0 .and. NumGlu(3).eq.0 .and. NumGlu(4).eq.0 ) return
!DEC$ IF (_DebugGeneralChecks==1)
      if( Quarks(2)%Mass.ne.0d0 .or. Quarks(3)%Mass.ne.0d0 ) print *, "CT insertions for massive quarks 2,3 not implemented"
!DEC$ ENDIF

      do n1a=0,NumGlu(1)
      do n3a=0,NumGlu(3)
         n1b = NumGlu(1)-n1a
         n3b = NumGlu(3)-n3a

         rIn =n1a+1
         rOut=NumGlu(1)+NumGlu(2)+n3a
         Eps2 = cur_g_2f(Gluons(rIn:rOut),Quarks(2:3),(/1+n1b+NumGlu(2)+n3a,n1b,NumGlu(2),n3a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         Eps2 = Eps2*PropFac1

         do n4a=0,NumGlu(4)
            n4b = NumGlu(4)-n4a
            ! Fer2
            rIn =NumGlu(1)+NumGlu(2)+n3a+1
            rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
            spiUDr(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(4:4),-Quarks(4)%PartType,(/n4a+n3b,n3b,n4a/))          ! undressed current
            spiDr(:) = (0d0,0d0)
            if(n3b.ge.1 .or. n4a.ge.1) then
               spiDr(:)  = cur_f_2f_massCT(Gluons(rIn:rOut),Quarks(4:4),-Quarks(4)%PartType,(/n4a+n3b,n3b,n4a/))   ! current dressed with mass CT's
               PMom2(:)  = Quarks(4)%Mom + SumMom(Gluons,rIn,rOut)
               PropFac2  = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(4)%Mass2)
               spiDr(:)  = ( spb2_(spiDr,PMom2)+Quarks(4)%Mass*spiDr(:) )*PropFac2
               spiDr(:)  = spiDr(:) + ( spb2_(spiUDr,PMom2)*(2d0*Quarks(4)%Mass) + (sc_(PMom2,PMom2)+Quarks(4)%Mass2)*spiUDr )*PropFac2**2
               spiUDr(:) = ( spb2_(spiUDr,PMom2)+Quarks(4)%Mass*spiUDr(:) )*PropFac2
               spiDr(:)  = vgq(eps2,spiDr)
            endif
            spiUDr(:) = vgq(eps2,spiUDr)

            PMom1 = Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom + SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a)
            if(n1a.ge.1 .or. n4b.ge.1) then
               PropFac1  = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(4)%Mass2)
               spiDr(:)  = ( spb2_(spiDr,PMom1)+Quarks(4)%Mass*spiDr(:) )*PropFac1
               spiDr(:)  = spiDr(:) + ( spb2_(spiUDr,PMom1)*(2d0*Quarks(4)%Mass) + (sc_(PMom1,PMom1)+Quarks(4)%Mass2)*spiUDr )*PropFac1**2
               spiUDr(:) = ( spb2_(spiUDr,PMom1)+Quarks(4)%Mass*spiUDr(:) )*PropFac1
            endif

            TmpExtRef = -1
            TmpQuark(1)%Mom  => PMom1(:)
            TmpQuark(1)%Pol  => spiDr(:)
            TmpQuark(1)%Mass => Quarks(4)%Mass
            TmpQuark(1)%Mass2=> Quarks(4)%Mass2
            TmpQuark(1)%ExtRef  => TmpExtRef
            TmpQuark(1)%PartType=> Quarks(4)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo

            tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),-TmpQuark(1)%PartType,(/counter-1,n1a,n4b/) )
            Res(:) = Res(:) + tmp(:)
            if( n1a.ge.1 .or. n4b.ge.1 ) then
              TmpQuark(1)%Pol  => spiUDr(:)
              tmp(:) = cur_f_2f_massCT(TmpGluons(1:counter-1),TmpQuark(1:1),-TmpQuark(1)%PartType,(/counter-1,n1a,n4b/) )
              Res(:) = Res(:) + tmp(:)
            endif
         enddo
      enddo
      enddo

   endif

return
END FUNCTION



































      recursive function g(e,k) result(res)
      implicit none
      complex(8), intent(in) :: e(:,:),k(:,:)
      complex(8)             :: e1(size(e,dim=1))
      complex(8)             :: res(size(e,dim=1))
      complex(8) :: k1(size(e,dim=1))
      complex(8) :: k2(size(e,dim=1))
      complex(8) :: k3(size(e,dim=1))
      complex(8) :: e2(size(e,dim=1))
      complex(8) :: e3(size(e,dim=1))
      complex(8) :: tmp(size(e,dim=1))
      complex(8) :: k1sq, k2sq, k3sq
      integer :: npart, m, m1

      npart = size(e,dim=2)

      if (npart == 1) then
         res = e(:,1)

      elseif (npart == 2) then
         res = vggg(e(:,1),k(:,1),e(:,2),k(:,2))
        else

           res = (0d0,0d0)

           do m=1,npart-1
              k1=sum(k(:,1:m),dim=2)
              k2=sum(k(:,m+1:npart),dim=2)

              e1 = g(e(:,1:m),k(:,1:m))
              e2 = g(e(:,m+1:npart),k(:,m+1:npart))

              tmp = vggg(e1,k1,e2,k2)

              if (m > 1) then
              k1sq = sc_(k1,k1)
               if (abs(k1sq) > propcut) then
              tmp = -(0d0,1d0)*tmp/k1sq
               else
                  tmp = (0d0,0d0)
               endif
              endif

              if (m + 1 < npart) then
              k2sq = sc_(k2,k2)
                if (abs(k2sq) > propcut) then
              tmp = -(0d0,1d0)*tmp/k2sq
                else
              tmp = (0d0,0d0)
                endif
              endif

              res = res + tmp

      if (m <= npart-2) then

        do m1=m+1,npart-1
!        e1 = g(e(:,1:m),k(:,1:m))  ! e1 is already computed
        e2=g(e(:,m+1:m1),k(:,m+1:m1))
        e3=g(e(:,m1+1:npart),k(:,m1+1:npart))
        k2 = sum(k(:,m+1:m1),dim=2)
        k3 = sum(k(:,m1+1:npart),dim=2)
        tmp = vgggg(e1,e2,e3)
        if (m > 1) then
        k1sq = sc_(k1,k1)
          if(abs(k1sq) > propcut) then
        tmp = -(0d0,1d0)*tmp/k1sq
         else
            tmp = (0d0,0d0)
          endif
        endif
        if (m+1 < m1) then
         k2sq = sc_(k2,k2)
           if(abs(k2sq) > propcut) then
           tmp = -(0d0,1d0)*tmp/k2sq
           else
            tmp = (0d0,0d0)
           endif
        endif
        if (m1+1 < npart) then
           k3sq = sc_(k3,k3)
           if(abs(k3sq) > propcut) then
           tmp = -(0d0,1d0)*tmp/k3sq
           else
           tmp = (0d0,0d0)
           endif
        endif
        res = res + tmp
        enddo

        endif
        enddo
        endif

      end function




      recursive function f(e,k,sp,p,mass,flout,flin,ms) result(res)
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
           sp2 = f(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,mass,flout,flin,ng1)
           if (ng1 >0.or.m>0) sp2 = spb2_(sp2,k2)+mass*sp2
           tmp = vqg(sp2,e1)

! print *, m,e1
! pause
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
           sp2=f(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,mass,flout,flin,ms1)

           if (ng2 > 0.or.m < ng1) sp2 = spb2_(sp2,k2)+mass*sp2
           tmp = vgq(e1,sp2)
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

      end function f



      recursive function bf(e,k,sp,p,mass,flout,flin,ms) result(res)
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
           sp2 = bf(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,mass,flout,flin,ng1)
           if (ng1>0.or.m>0) sp2 = spi2_(k2,sp2)+mass*sp2

            tmp = vbqg(sp2,e1)

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
           sp2=bf(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,mass,flout,flin,ms1)

           if (ng2 > 0.or.m < ng1) sp2 = spi2_(k2,sp2)+mass*sp2

           tmp = vgbq(e1,sp2)

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
      end function bf






      recursive function fmCT(e,k,sp,p,mass,flout,flin,ms,CTIns) result(res)
      implicit none
      logical CTIns
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
           sp2 = fmCT(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,mass,flout,flin,ng1,.false.)

           if (ng1 >0.or.m>0) then
               if (CTIns) then
                  sp2 = spb2_(sp2,k2)*(2d0*mass)+ (sc_(k2,k2)+mass**2)*sp2

                  if (abs(k2sq) > propcut) then
                        sp2 =  sp2/k2sq * ((0d0,1d0))
                  else
                        sp2 = (0d0,0d0)
                  endif

                  if(m>1 .or. ng1>0) then
                     sp3 = fmCT(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,mass,flout,flin,ng1,.true.)
                     sp3 = spb2_(sp3,k2)+mass*sp3
                     sp2(:) = sp2(:) + sp3(:)
                  endif
               else
                  sp2 = spb2_(sp2,k2)+mass*sp2
               endif
           else
               if(CTIns) cycle
           endif

           tmp = vqg(sp2,e1)

           if (ng1>0.or.m>0) then
              if (abs(k2sq) > propcut) then
                  tmp =  (0d0,1d0)/k2sq*tmp
              else
                  tmp = (0d0,0d0)
               endif
           endif

           if (m < ng2-1)  then
              if(abs(k1sq) > propcut) then
                  tmp = -(0d0,1d0)/k1sq*tmp
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
           sp2=fmCT(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,mass,flout,flin,ms1,.false.)

           if (ng2 > 0.or.m < ng1) then
               if (CTIns) then
                  sp2 = spb2_(sp2,k2)*(2d0*mass)+ (sc_(k2,k2)+mass**2)*sp2

                  if (abs(k2sq) > propcut) then
                        sp2 =  sp2/k2sq * ((0d0,1d0))
                  else
                        sp2 = (0d0,0d0)
                  endif

                  if(m<ng1-1 .or. ng2>0) then
                     sp3 = fmCT(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,mass,flout,flin,ms1,.true.)
                     sp3 = spb2_(sp3,k2)+mass*sp3
                     sp2(:) = sp2(:) + sp3(:)
                  endif
               else
                  sp2 = spb2_(sp2,k2)+mass*sp2
               endif
           else
               if(CTIns) cycle
           endif

           tmp = vgq(e1,sp2)

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


         return
      end function fmCT







! FUNCTION cur_f_4f_massCT_old(Gluons,Quarks,Quark1PartType,NumGlu) result(res)
! implicit none
! integer :: NumGlu(0:4),Quark1PartType
! type(PtrToParticle) :: Gluons(1:),Quarks(2:4)            ! off-shell quark is not included in Quarks(.)
! integer,target :: TmpExtRef,TmpPartType
! complex(8) :: res(1:Ds),tmp(1:Ds)
! complex(8),target :: spiDr(1:Ds),spiUDr(1:Ds)
! complex(8) :: eps1(1:Dv)
! complex(8) :: eps2(1:Dv)
! type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(4)),TmpQuarkDr(1:1),TmpQuarkUDr(1:1)
! complex(8) :: PropFac1,PropFac2
! complex(8),target :: pmom1(1:Dv)
! complex(8) :: pmom2(1:Dv)
! integer :: n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b
! integer :: rIn,rOut,i,counter
!
!
!
! !DEC$ IF (_DebugGeneralChecks==1)
!     if( Quarks(3)%Mass.ne.0d0 .or. Quarks(4)%Mass.ne.0d0 ) print *, "CT insertions for massive quarks 3,4 not implemented"
! !DEC$ ENDIF
!
!    Res(:)=(0d0,0d0)
!    if( NumGlu(1).eq.0 .and. NumGlu(2).eq.0 .and. NumGlu(4).eq.0 ) return
!
!    if( NumGlu(1).eq.1 ) then
!          Eps2 = cur_g_2f(Gluons(1:0),Quarks(3:4),(/1,0,0,0/))
!          PMom1(:) = Quarks(3)%Mom + Quarks(4)%Mom
!          PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
!          Eps2 = Eps2*PropFac1
!
!             spiUDr(:) = cur_f_2f(Gluons(1:1),Quarks(2:2),-Quarks(2)%PartType,(/1,1,0/))
!             PMom2(:) = Quarks(2)%Mom + SumMom(Gluons,1,1)
!             PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(2)%Mass2)
!             spiDr(:)  = ( spb2_(spiUDr,PMom2)*(2d0*Quarks(2)%Mass) + (sc_(PMom2,PMom2)+Quarks(2)%Mass2)*spiUDr )*PropFac2**2
!             spiDr(:)  = vqg(spiDr, eps2)
!             Res(:) = Res(:) + spiDr(:)
!
!             spiUDr(:) = vqg(Quarks(2)%Pol(:),eps2)
!             PMom1 = Quarks(2)%Mom+Quarks(3)%Mom+Quarks(4)%Mom
!             PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(2)%Mass2)
!             spiDr(:)  = ( spb2_(spiUDr,PMom1)*(2d0*Quarks(2)%Mass) + (sc_(PMom1,PMom1)+Quarks(2)%Mass2)*spiUDr )*PropFac1**2
!             spiDr(:) = vgq(Gluons(1)%Pol,spiDr)
!             Res(:) = Res(:) + spiDr(:)
!
!    elseif( NumGlu(2).eq.1 ) then
!
!          Eps2 = cur_g_2f(Gluons(1:0),Quarks(3:4),(/1,0,0,0/))
!          PMom1(:) = Quarks(3)%Mom + Quarks(4)%Mom
!          PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
!          Eps2 = Eps2*PropFac1
!
!             spiUDr(:) = cur_f_2f(Gluons(1:1),Quarks(2:2),-Quarks(2)%PartType,(/1,0,1/))
!             PMom2(:) = Quarks(2)%Mom + SumMom(Gluons,1,1)
!             PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(2)%Mass2)
!             spiDr(:)  = ( spb2_(spiUDr,PMom2)*(2d0*Quarks(2)%Mass) + (sc_(PMom2,PMom2)+Quarks(2)%Mass2)*spiUDr )*PropFac2**2
!             spiDr(:)  = vqg(spiDr, eps2)
!             Res(:) = Res(:) + spiDr(:)
!     endif
! return
! END FUNCTION









! ---------------------------- currents with a W boson coupling to a top-bot-line -----------------------------------------
!------------------------------           works only in 4D because of Chir        -----------------------------------------


FUNCTION cur_f_2fW(Gluons,Quarks,Boson,NumGlu) result(Res)           ! Quarks(:) DOES include the OFF-shell quark, in contrast to all other routines!
implicit none
complex(8) :: Res(1:Ds)
integer :: NumGlu(0:2),i,rIn,rOut,Quark1PartType
type(PtrToParticle) :: Gluons(1:),Boson,Quarks(1:2) 
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
      call Error("This Flavor is not allowed in cur_f_2fW",Quark1PartType)
   endif

   if( abs(Quarks(2)%PartType).eq.4 ) then!  4=Str_
     FerFla2="str"
   else
      call Error("This Flavor is not allowed in cur_f_2fW",Quarks(2)%PartType)
   endif

   if( NumGlu(0)-NumGlu(1)-NumGlu(2).ne.0 ) call Error("Wrong NumGlu in cur_f_2fV",NumGlu(0)-NumGlu(1)-NumGlu(2))

   rIn =1
   rOut=NumGlu(0)
   if( Quarks(2)%PartType .gt.0 ) then      !    X----->----
      Res(:) = fW(GluPol(1:Dv,rIn:rOut),GluMom(1:Dv,rIn:rOut),QuarkPol(1:Ds),QuarkMom(1:Dv),Quarks(1)%Mass,FerFla2,FerFla1,Boson%Pol(1:Dv),Boson%Mom(1:Dv),NumGlu(1))
   else                                     !    X-----<----
      Res(:) = bfW(GluPol(1:Dv,rIn:rOut),GluMom(1:Dv,rIn:rOut),QuarkPol(1:Ds),QuarkMom(1:Dv),Quarks(1)%Mass,FerFla2,FerFla1,Boson%Pol(1:Dv),Boson%Mom(1:Dv),NumGlu(1))
   endif


return
END FUNCTION





      recursive function fW(e,k,sp,p,mass,flout,flin,eW,kW,ms) result(res)
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


      if (ng2 < 0) write(*,*) 'WRONG DEFINITION OF CURRENT fW'

if (ngluon == 0) then
         res = vbqW(sp,eW)
else

       res = (0d0,0d0)
       do m=0,ng2-1
           k1 = sum(k(:,ng1+1+m:ngluon),dim=2)
           e1=g(e(:,ng1+1+m:ngluon),k(:,ng1+1+m:ngluon))
           k1sq=sc_(k1,k1)

           k2 = sum(k(:,1:ng1+m),dim=2)
           k2 = k2 + p + kW
           k2sq = sc_(k2,k2)-mass**2
           sp2 = fW(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,mass,flout,flin,eW,kW,ng1)
           sp2 = spb2_(sp2,k2)+mass*sp2

           tmp = vqg(sp2,e1)
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
           sp2=fW(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,mass,flout,flin,eW,kW,ms1)

           sp2 = spb2_(sp2,k2)+mass*sp2

           tmp = vgq(e1,sp2)
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



        sp2 = f(e,k,sp,p,0d0,flout,flout,ms) !  mass set to zero because this is a bottom quark
        k2 = sum(k(:,1:ngluon),dim=2)
        k2 = k2 + p
        k2sq = sc_(k2,k2)   !-mass**2

        sp2 = spb2_(sp2,k2) !+mass*sp2

        tmp = vbqW(sp2,eW)
        if (abs(k2sq) > propcut) then
               tmp = (0d0,1d0)/k2sq*tmp
        else
                tmp = (0d0,0d0)
        endif
        res = res + tmp

endif

end function fW




      recursive function bfW(e,k,sp,p,mass,flout,flin,eW,kW,ms) result(res)
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
         res = vWq(eW,sp)
else

       res = (0d0,0d0)
       do m=0,ng2-1
           k1 = sum(k(:,ng1+1+m:ngluon),dim=2)
           e1=g(e(:,ng1+1+m:ngluon),k(:,ng1+1+m:ngluon))
           k1sq=sc_(k1,k1)

           k2 = sum(k(:,1:ng1+m),dim=2)
           k2 = -k2 - p - kW
           k2sq = sc_(k2,k2)-mass**2

           sp2 = bfW(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,mass,flout,flin,eW,kW,ng1)
           sp2 = spi2_(k2,sp2)+mass*sp2

           tmp = vbqg(sp2,e1)

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
           sp2=bfW(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,mass,flout,flin,eW,kW,ms1)

           sp2 = spi2_(k2,sp2)+mass*sp2

           tmp = vgbq(e1,sp2)

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



        sp2 = bf(e,k,sp,p,0d0,flout,flout,ms) !  mass set to zero because this is a bottom quark
        k2 = sum(k(:,1:ngluon),dim=2)
        k2 = -k2 - p
        k2sq = sc_(k2,k2)   !-mass**2

        sp2 = spb2_(sp2,k2) !+mass*sp2

        tmp = vWq(eW,sp2)
        if (abs(k2sq) > propcut) then
               tmp = (0d0,1d0)/k2sq*tmp
        else
               tmp = (0d0,0d0)
        endif
        res = res + tmp

endif

      end function bfW





! ---------------------------- currents with a vector boson coupling to a quark-line -----------------------------------------
!------------------------------           works only in 4D because of Chir        -----------------------------------------


FUNCTION cur_f_2fV(Gluons,Quarks,Boson,NumGlu) result(Res)           ! Quarks(:) DOES include the OFF-shell quark, in contrast to all other routines!
implicit none
complex(8) :: Res(1:Ds)
integer :: NumGlu(0:2),i,rIn,rOut,Quark1PartType
type(PtrToParticle) :: Gluons(1:),Boson,Quarks(1:2) 
complex(8) :: GluMom(1:Dv,NumGlu(0)), QuarkMom(1:Dv)
complex(8) :: GluPol(1:Dv,NumGlu(0)), QuarkPol(1:Ds)



   do i=1,NumGlu(0)
    GluMom(1:Dv,i) = Gluons(i)%Mom(1:Dv)
    GluPol(1:Dv,i) = Gluons(i)%Pol(1:Dv)
   enddo
   QuarkMom(1:Dv) = Quarks(2)%Mom(1:Dv)
   QuarkPol(1:Ds) = Quarks(2)%Pol(1:Ds)

   if( Quarks(1)%PartType.ne.-Quarks(2)%PartType ) call Error("Wrong quark flavors in cur_f_2fV")
   if( NumGlu(0)-NumGlu(1)-NumGlu(2).ne.0 ) call Error("Wrong NumGlu in cur_f_2fV",NumGlu(0)-NumGlu(1)-NumGlu(2))

   rIn =1
   rOut=NumGlu(0)
   if( Quarks(2)%PartType .gt.0 ) then      !    X----->----
      Res(:) = fV(GluPol(1:Dv,rIn:rOut),GluMom(1:Dv,rIn:rOut),QuarkPol(1:Ds),QuarkMom(1:Dv),Quarks(1)%Mass,Quarks(1)%PartType,Boson%Pol(1:Dv),Boson%Mom(1:Dv),NumGlu(1))
   else                                     !    X-----<----
      Res(:) = bfV(GluPol(1:Dv,rIn:rOut),GluMom(1:Dv,rIn:rOut),QuarkPol(1:Ds),QuarkMom(1:Dv),Quarks(1)%Mass,Quarks(1)%PartType,Boson%Pol(1:Dv),Boson%Mom(1:Dv),NumGlu(1))
   endif


return
END FUNCTION







      recursive function fV(e,k,sp,p,mass,QuarkFlavor,eV,kV,ms) result(res)
      use ModParameters
      implicit none
      complex(8), intent(in) :: e(:,:), k(:,:)
      complex(8), intent(in) :: sp(:), p(:)
      complex(8), intent(in) :: eV(:), kV(:)
      integer, intent(in) ::  ms,QuarkFlavor
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
      real(8) :: mass,couplVQQ_left,couplVQQ_right
      character,parameter :: FerFla*3="dum" ! dummy, only used for check of flavor consistency inside the functions f,bf


      ngluon = size(e,dim=2)
      ng1 = ms   !#gluons to the left of a f-line
      ng2 = ngluon - ms  !#gluons to the right of the f-line

      if (ng2 < 0) write(*,*) 'WRONG DEFINITION OF CURRENT fV'

      if( abs(QuarkFlavor).eq.Up_ .or. abs(QuarkFlavor).eq.Chm_ ) then 
          couplVQQ_left  = couplZUU_left
          couplVQQ_right = couplZUU_right
      elseif( abs(QuarkFlavor).eq.Dn_ .or. abs(QuarkFlavor).eq.Str_  ) then 
          couplVQQ_left  = couplZUU_left
          couplVQQ_right = couplZUU_right
      elseif( abs(QuarkFlavor).eq.Top_ .or. abs(QuarkFlavor).eq.Bot_) then!   note that Bot_ is treated as top quark in TOPAZ!
          couplVQQ_left  = couplZTT_left
          couplVQQ_right = couplZTT_right
      endif
if (ngluon == 0) then
         res = vbqV(sp,eV,couplVQQ_left,couplVQQ_right)
else

       res = (0d0,0d0)
       do m=0,ng2-1
           k1 = sum(k(:,ng1+1+m:ngluon),dim=2)
           e1=g(e(:,ng1+1+m:ngluon),k(:,ng1+1+m:ngluon))
           k1sq=sc_(k1,k1)

           k2 = sum(k(:,1:ng1+m),dim=2)
           k2 = k2 + p + kV
           k2sq = sc_(k2,k2)-mass**2
           sp2 = fV(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,mass,QuarkFlavor,eV,kV,ng1)
           sp2 = spb2_(sp2,k2)+mass*sp2

           tmp = vqg(sp2,e1)
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
           k2 = k2 + p + kV
           k2sq = sc_(k2,k2) - mass**2
           ms1 = ng1 - m
           sp2=fV(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,mass,QuarkFlavor,eV,kV,ms1)

           sp2 = spb2_(sp2,k2)+mass*sp2

           tmp = vgq(e1,sp2)
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


        sp2 = f(e,k,sp,p,mass,FerFla,FerFla,ms)
        k2 = sum(k(:,1:ngluon),dim=2)
        k2 = k2 + p
        k2sq = sc_(k2,k2)  - mass**2

        sp2 = spb2_(sp2,k2)+ mass*sp2

        tmp = vbqV(sp2,eV,couplVQQ_left,couplVQQ_right)
        if (abs(k2sq) > propcut) then
               tmp = (0d0,1d0)/k2sq*tmp
        else
                tmp = (0d0,0d0)
        endif
        res = res + tmp

endif

end function fV




      recursive function bfV(e,k,sp,p,mass,QuarkFlavor,eV,kV,ms) result(res)
      use ModParameters
      implicit none
      complex(8), intent(in) :: e(:,:), k(:,:)
      complex(8), intent(in) :: sp(:), p(:)
      complex(8), intent(in) :: eV(:), kV(:)
      integer, intent(in) ::  ms,QuarkFlavor
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
      real(8) :: mass,couplVQQ_left,couplVQQ_right
      character,parameter :: FerFla*3="dum" ! dummy, only used for check of flavor consistency inside the functions f,bf


      ngluon = size(e,dim=2)
      ng1 = ms   !#gluons to the left of a f-line
      ng2 = ngluon - ms  !#gluons to the right of the f-line

      if (ng2 < 0) write(*,*) 'WRONG DEFINITION OF CURRENT fbV'

      if( abs(QuarkFlavor).eq.Up_ .or. abs(QuarkFlavor).eq.Chm_ ) then 
          couplVQQ_left  = couplZUU_left
          couplVQQ_right = couplZUU_right
      elseif( abs(QuarkFlavor).eq.Dn_ .or. abs(QuarkFlavor).eq.Str_  ) then 
          couplVQQ_left  = couplZUU_left
          couplVQQ_right = couplZUU_right
      elseif( abs(QuarkFlavor).eq.Top_ .or. abs(QuarkFlavor).eq.Bot_) then!   note that Bot_ is treated as top quark in TOPAZ!
          couplVQQ_left  = couplZTT_left
          couplVQQ_right = couplZTT_right
      endif

if (ngluon == 0) then
         res = vVq(eV,sp,couplVQQ_left,couplVQQ_right)
else

       res = (0d0,0d0)
       do m=0,ng2-1
           k1 = sum(k(:,ng1+1+m:ngluon),dim=2)
           e1=g(e(:,ng1+1+m:ngluon),k(:,ng1+1+m:ngluon))
           k1sq=sc_(k1,k1)

           k2 = sum(k(:,1:ng1+m),dim=2)
           k2 = -k2 - p - kV
           k2sq = sc_(k2,k2)-mass**2

           sp2 = bfV(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,mass,QuarkFlavor,eV,kV,ng1)
           sp2 = spi2_(k2,sp2)+mass*sp2

           tmp = vbqg(sp2,e1)

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
           k2 = -k2 - p - kV
           k2sq = sc_(k2,k2) - mass**2
           ms1 = ng1 - m
           sp2=bfV(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,mass,QuarkFlavor,eV,kV,ms1)

           sp2 = spi2_(k2,sp2)+mass*sp2

           tmp = vgbq(e1,sp2)

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



        sp2 = bf(e,k,sp,p,mass,FerFla,FerFla,ms)
        k2 = sum(k(:,1:ngluon),dim=2)
        k2 = -k2 - p
        k2sq = sc_(k2,k2)   -mass**2

        sp2 = spb2_(sp2,k2) +mass*sp2

        tmp = vVq(eV,sp2,couplVQQ_left,couplVQQ_right)
        if (abs(k2sq) > propcut) then
               tmp = (0d0,1d0)/k2sq*tmp
        else 
               tmp = (0d0,0d0)
        endif
        res = res + tmp

endif

      end function bfV








! ---------------------------- currents with scalars coupling to gluons  -----------------------------------------




!     assuming the color flow:   s = x--->----
      RECURSIVE FUNCTION cur_s_2s(Gluons,Scalar,NumGlu) result(res)! gauge invariance checked for up to 4 gluons
      implicit none
      type(PtrToParticle) :: Gluons(1:),Scalar(2:2)      ! off-shell scalar is not included
      integer ::  NumGlu(0:2)
      integer :: ms1,m1,m2,ng1, ng2, ngluon
      complex(8) :: res,tmp,res1,res2,res3,res4,res5
      complex(8) :: sc2
      complex(8) :: k1(1:Dv),k2(1:Dv),k3(1:Dv)
      complex(8) :: e1(1:Dv)
      complex(8) :: e2(1:Dv)
      complex(8) :: k1sq,k2sq,k3sq



       ngluon = NumGlu(0)
       ng1 = NumGlu(1)           !#gluons to the left of a s-line
       ng2 = NumGlu(2)           !#gluons to the right of the s-line

!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2).ne.0 ) print *, "wrong number of gluons in cur_s_2s"
!DEC$ ENDIF

       if (ngluon .eq. 0) then
         res = Scalar(2)%Pol(1)
         return
       endif
       res = (0d0,0d0)
res1=(0d0,0d0); res2=(0d0,0d0); res3=(0d0,0d0); res4=(0d0,0d0); res5=(0d0,0d0)


!    s-s-g current to the right
       do m1=0,ng2-1
           k1 = SumMom(Gluons,ng1+1+m1,ngluon)
           e1 = cur_g(Gluons(ng1+1+m1:ngluon),ngluon-m1-ng1+1)
           k1sq=sc_(k1,k1)

           k2 = SumMom(Gluons,1,ng1+m1)
           k2 = k2 + Scalar(2)%Mom
           k2sq = sc_(k2,k2) - Scalar(2)%Mass2
           sc2 = cur_s_2s(Gluons(1:ng1+m1),Scalar(2:2),(/ng1+m1,ng1,m1/))
           if( Scalar(2)%PartType.gt.0 ) then
              tmp = csg(e1,k1,k2) * sc2
           else
              tmp = cbsg(e1,k1,k2) * sc2
           endif

           if (m1 < ng2-1)  then
              if(abs(k1sq) > propcut) then
                  tmp = -(0d0,1d0)/k1sq*tmp! gluon prop
              else
                  tmp = (0d0,0d0)
              endif
           endif
           if (ng1>0.or.m1>0) then
              if (abs(k2sq) > propcut) then
                  tmp =  (0d0,1d0)/k2sq*tmp! scalar prop
              else
                  tmp = (0d0,0d0)
               endif
           endif
           res = res + tmp
res1=res1+tmp
        enddo




!     s-g-s current to the left
        do m1=1,ng1
!  print *, "loop 2",ng1,ng2
           k1 = SumMom(Gluons,1,m1)
           e1 = cur_g(Gluons(1:m1),m1+1)
           k1sq = sc_(k1,k1)

           k2 = SumMom(Gluons,m1+1,ngluon)
           k2 = k2 + Scalar(2)%Mom
           k2sq = sc_(k2,k2) - Scalar(2)%Mass2
           ms1 = ng1 - m1
           sc2 = cur_s_2s(Gluons(m1+1:ngluon),Scalar(2:2),(/ngluon-m1,ms1,ngluon-ng1/))
           if( Scalar(2)%PartType.gt.0 ) then
              tmp = cgs(e1,k1,k2) * sc2
           else
              tmp = cgbs(e1,k1,k2) * sc2
           endif

           if (m1 > 1) then
              if (abs(k1sq) > propcut) then
                  tmp=-(0d0,1d0)/k1sq*tmp
              else
                  tmp = (0d0,0d0)
              endif
           endif

           if (ng2 > 0.or. m1<ng1) then
            if (abs(k2sq) > propcut) then
              tmp=(0d0,1d0)/k2sq*tmp
            else
              tmp = (0d0,0d0)
            endif
           endif
           res = res + tmp
res2=res2+tmp
           enddo




!     s-s-g-g current to the right
        do m2=2,ng2
! print *, "loop 3",ng1,ng2
           k3 = SumMom(Gluons,1,ngluon-m2)
           k3 = k3 + Scalar(2)%Mom
           k3sq = sc_(k3,k3) - Scalar(2)%Mass2
           sc2 = cur_s_2s(Gluons(1:ngluon-m2),Scalar(2:2),(/ngluon-m2,ng1,ngluon-m2-ng1/))
           if( ngluon-m2-ng1>0 .or. ng1>0) then! here was a bug! here was another bug
             if (abs(k3sq) > propcut) then
               sc2= (0d0,1d0)/k3sq*sc2
             else
               sc2 = (0d0,0d0)
               cycle
             endif
           endif

           do m1=1,m2-1
              k1 = SumMom(Gluons,ngluon-m2+1,ngluon-m2+m1)
              e1 = cur_g(Gluons(ngluon-m2+1:ngluon-m2+m1),m1+1)
              k1sq = sc_(k1,k1)

              k2 = SumMom(Gluons,ngluon-m2+m1+1,ngluon)
              e2 = cur_g(Gluons(ngluon-m2+m1+1:ngluon),m2-m1+1)
              k2sq = sc_(k2,k2)
              tmp = csgg(e1,e2) * sc2
              if( m1>1 ) then
                if (abs(k1sq) > propcut) then
                   tmp=-(0d0,1d0)/k1sq*tmp
                else
                   tmp = (0d0,0d0)
                endif
              endif
              if( m2-m1>1 ) then
                if (abs(k2sq) > propcut) then
                   tmp=-(0d0,1d0)/k2sq*tmp
                else
                   tmp = (0d0,0d0)
                endif                                      
              endif
              res = res + tmp
res3=res3+tmp
           enddo
        enddo



!     s-g-g-s current to the left
        do m2=2,ng1
! print *, "loop 4",ng1,ng2
           k3 = SumMom(Gluons,m2+1,ngluon)
           k3 = k3 + Scalar(2)%Mom
           k3sq = sc_(k3,k3) - Scalar(2)%Mass2
           ms1 = ng1 - m2
           sc2 = cur_s_2s(Gluons(m2+1:ngluon),Scalar(2:2),(/ngluon-m2,ms1,ngluon-ng1/))
           if( ms1>0 .or. ngluon-ng1>0) then! here was a bug
           if (abs(k3sq) > propcut) then
               sc2= (0d0,1d0)/k3sq*sc2
           else
               sc2 = (0d0,0d0)
               cycle
           endif
           endif
           do m1=1,m2-1
              k1 = SumMom(Gluons,1,m1)
              e1 = cur_g(Gluons(1:m1),m1+1)
              k1sq = sc_(k1,k1)

              k2 = SumMom(Gluons,m1+1,m2)
              e2 = cur_g(Gluons(m1+1:m2),m2-m1+1)
              k2sq = sc_(k2,k2)

              tmp = cggs(e1,e2) * sc2
              if( m1>1 ) then
                if (abs(k1sq) > propcut) then
                   tmp=-(0d0,1d0)/k1sq*tmp
                else
                   tmp = (0d0,0d0)
                endif                   
              endif
              if( m2-m1>1 ) then
                if (abs(k2sq) > propcut) then
                   tmp=-(0d0,1d0)/k2sq*tmp
                else
                   tmp = (0d0,0d0)
                endif                                      
              endif
              res = res + tmp
res4=res4+tmp
           enddo
        enddo



!     s-g-s-g current to the left and right
      if( ng1.ge.1 .and. ng2.ge.1 ) then
! print *, "loop 5",ng1,ng2
        do m1=1,ng1
           k1 = SumMom(Gluons,1,m1)
           e1 = cur_g(Gluons(1:m1),m1+1)
           if( m1.gt.1 ) then 
              k1sq = sc_(k1,k1)
              if (abs(k1sq) > propcut) then
                  e1 = -(0d0,1d0)/k1sq * e1
              else
                  e1 = (0d0,0d0)
                  cycle
              endif
           endif
           do m2=0,ng2-1
! print *, "loop 5 ", ng1,ng2
! print *, "loop 5 ", m1,m2
              k2 = SumMom(Gluons,ng1+m2+1,ngluon)
              e2 = cur_g(Gluons(ng1+m2+1:ngluon),ngluon-ng1-m2-1+1+1)
              if( ngluon-ng1-m2-1+1.gt.1 ) then
                  k2sq = sc_(k2,k2)
                  if (abs(k2sq) > propcut) then
                      e2 = -(0d0,1d0)/k2sq * e2
                  else
                      e2 = (0d0,0d0)
                      cycle
                  endif
              endif 
              sc2 = cur_s_2s(Gluons(m1+1:ng1+m2),Scalar(2:2),(/ng1-m1+m2,ng1-m1,m2/))
              if( ng1-m1.gt.0 .or. m2.gt.0 ) then
                k3 = SumMom(Gluons,m1+1,ng1+m2)
                k3 = k3 + Scalar(2)%Mom
                k3sq = sc_(k3,k3) - Scalar(2)%Mass2
                if (abs(k3sq) > propcut) then
                   sc2 = (0d0,1d0)/k3sq*sc2
                else
                   sc2 = (0d0,0d0)
                endif
              endif
              tmp = cgsg(e1,e2) * sc2
              res = res + tmp
res5=res5+tmp
           enddo
        enddo
      endif

! if(ng1.eq.2 .and. ng2.eq.2) then
!  print *,"s_2s", ng1,ng2
!  print *,res1
!  print *,res2
!  print *,res3
!  print *,res4
!  print *,res5
! pause
! endif
! print *, "exit",ngluon,ng1,ng2
      end function cur_s_2s




FUNCTION cur_s_ssff(Gluons,Scalar,Quarks,NumGlu) result(res)!  checked gauge invariance for 1 gluon
implicit none
integer :: NumGlu(0:4)
type(PtrToParticle) :: Gluons(1:),Quarks(3:4),Scalar(2:2)      ! off-shell scalar is not included
integer,target :: TmpExtRef
complex(8) :: res,tmp
complex(8),target :: Sca0(1:1)
complex(8) :: Sca1
complex(8) :: eps1(1:Dv),eps2(1:Dv),EpsX(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(4)),TmpScalar(1:1)
complex(8) :: PropFac1,PropFac2,PROPFAC4
complex(8),target :: pmom1(1:Dv)
complex(8) :: pmom2(1:Dv),pmom4(1:Dv)
integer :: n1a,n1b,n1c,n2a,n2b,n2c,n3a,n3b,n4a,n4b,n4c
integer :: rIn,rOut,i,counter


!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4).ne.0 ) print *, "wrong number of gluons in cur_s_ssff"
!DEC$ ENDIF

      Res=(0d0,0d0)
      do n1a=0,NumGlu(1)
      do n2a=0,NumGlu(2)
      do n4a=0,NumGlu(4)
      do n1c=0,NumGlu(1)-n1a
      do n2c=0,NumGlu(2)-n2a
      do n4c=0,NumGlu(4)-n4a

         n1b = NumGlu(1)-n1a-n1c
         n2b = NumGlu(2)-n2a-n2c
         n4b = NumGlu(4)-n4a-n4c
         if( n1c.gt.0 .and. (n2c+n4c).gt.0  ) cycle
         if( n2c.gt.0 .and. (n1c+n4c).gt.0  ) cycle
         if( n4c.gt.0 .and. (n1c+n2c).gt.0  ) cycle

!print *, "n1",n1a,n1b,n1c
!print *, "n2",n2a,n2b,n2c
!print *, "n4",n4a,n4b,n4c

         rIn =NumGlu(1)+n2a+n2c+1
         rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
         Eps2 = cur_g_2f(Gluons(rIn:rOut),Quarks(3:4),(/1+n2b+NumGlu(3)+n4a,n2b,NumGlu(3),n4a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom + Quarks(4)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
         Eps2 = Eps2*PropFac1

            rIn =n1a+n1c+1
            rOut=NumGlu(1)+n2a
            Sca1 = cur_s_2s(Gluons(rIn:rOut),Scalar(2:2),(/n2a+n1b,n1b,n2a/))
!if(n1b.eq.1) print *, "xx",(0d0,1d0)/dsqrt(2d0)*2d0*(Scalar(2)%Mom.dot.Gluons(1)%Pol) & 
!                          *(0d0,1d0)/2d0/(Scalar(2)%Mom.dot.Gluons(1)%Mom) & 
!                          *(0d0,-1d0)/dsqrt(2d0)*(2d0*(Scalar(2)%Mom+Gluons(1)%Mom).dot.Eps2) 
            PMom2(:) = Scalar(2)%Mom + SumMom(Gluons,rIn,rOut)
            if(n1b.ge.1 .or. n2a.ge.1) then
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Scalar(2)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Scalar(2)%Mass2).lt.PropCut ) cycle
               Sca1 = Sca1*PropFac2
            endif

            if( n1c.gt.0 ) then
                rIn =n1a+1
                rOut=n1a+n1c
                EpsX(:) = cur_g(Gluons(rIn:rOut),1+n1c)
                if(n1c.gt.1) then
                    PMom4 = SumMom(Gluons,rIn,rOut)
                    PropFac4 = (0d0,-1d0)/sc_(PMom4,PMom4)
                    if( abs(sc_(PMom4,PMom4)).lt.PropCut ) cycle
                    EpsX(:) = EpsX(:)*PropFac4
                endif
                Sca0 = cgsg(EpsX,Eps2) * Sca1!   here was a bug
            elseif( n2c.gt.0) then
                rIn =NumGlu(1)+n2a+1
                rOut=NumGlu(1)+n2a+n2c
                EpsX(:) = cur_g(Gluons(rIn:rOut),1+n2c)
                if(n2c.gt.1) then
                    PMom4 = SumMom(Gluons,rIn,rOut)
                    PropFac4 = (0d0,-1d0)/sc_(PMom4,PMom4)
                    if( abs(sc_(PMom4,PMom4)).lt.PropCut ) cycle
                    EpsX(:) = EpsX(:)*PropFac4
                endif
                Sca0 = csgg(Eps2,EpsX) * Sca1
            elseif( n4c.gt.0 ) then
                rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+n4c
                EpsX(:) = cur_g(Gluons(rIn:rOut),1+n4c)
                if(n4c.gt.1) then
                    PMom4 = SumMom(Gluons,rIn,rOut)
                    PropFac4 = (0d0,-1d0)/sc_(PMom4,PMom4)
                    if( abs(sc_(PMom4,PMom4)).lt.PropCut ) cycle
                    EpsX(:) = EpsX(:)*PropFac4
                endif
                Sca0 = csgg(Eps2,EpsX) * Sca1
            else
                  if( Scalar(2)%PartType.gt.0 ) then
                      Sca0 = csg(Eps2,PMom1,PMom2) * Sca1
                  else
                      Sca0 = cbsg(Eps2,PMom1,PMom2) * Sca1
                  endif
            endif


            PMom1 = Scalar(2)%Mom+Quarks(3)%Mom+Quarks(4)%Mom+SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+n4c)
            if(n1a.ge.1 .or. n4b.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Scalar(2)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Scalar(2)%Mass2).lt.PropCut ) cycle
               Sca0 = Sca0*PropFac1
            endif


            TmpScalar(1)%Mom  => PMom1(:)
            TmpScalar(1)%Pol  => Sca0
            TmpScalar(1)%Mass => Scalar(2)%Mass
            TmpScalar(1)%Mass2=> Scalar(2)%Mass2
            TmpExtRef = -1
            TmpScalar(1)%ExtRef => TmpExtRef
            TmpScalar(1)%PartType => Scalar(2)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+n4c+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp = cur_s_2s(TmpGluons(1:counter-1),TmpScalar(1:1),(/counter-1,n1a,n4b/) )
!if(n1a.eq.1) print *, "yy",(0d0,-1d0)/dsqrt(2d0)*2d0*(Scalar(2)%Mom.dot.Eps2) & 
!                          *(0d0,1d0)/2d0/((Scalar(2)%Mom.dot.Quarks(3)%Mom)+(Scalar(2)%Mom.dot.Quarks(4)%Mom)+(Quarks(4)%Mom.dot.Quarks(3)%Mom)) &
!                          *(0d0,1d0)/dsqrt(2d0)*2d0*((Scalar(2)%Mom+Quarks(3)%Mom+Quarks(4)%Mom).dot.Gluons(1)%Pol)


! print *, "before",res
            Res = Res + tmp
! print *,"add",(tmp)
! print *,"sum",(res)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!print *,"sum",(res)


return
END FUNCTION





FUNCTION cur_s_sffs(Gluons,Scalar,Quarks,NumGlu) result(res)!  checked gauge invariance for 1 gluon
implicit none
integer :: NumGlu(0:4)
type(PtrToParticle) :: Gluons(1:),Quarks(2:3),Scalar(4:4)      ! off-shell scalar is not included
integer,target :: TmpExtRef
complex(8) :: res,tmp
complex(8),target :: Sca0(1:1)
complex(8) :: Sca1
complex(8) :: eps1(1:Dv),eps2(1:Dv),epsX(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(4)),TmpScalar(1:1)
complex(8) :: PropFac1,PropFac2,PropFac4
complex(8),target :: pmom1(1:Dv)
complex(8) :: pmom2(1:Dv),pmom4(1:Dv)
integer :: n1a,n1b,n1c,n3a,n3b,n3c,n4a,n4b,n4c
integer :: rIn,rOut,i,counter


!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4).ne.0 ) print *, "wrong number of gluons in cur_s_sffs"
!DEC$ ENDIF

      Res=(0d0,0d0)
      do n1a=0,NumGlu(1)
      do n3a=0,NumGlu(3)
      do n4a=0,NumGlu(4)
      do n1c=0,NumGlu(1)-n1a
      do n3c=0,NumGlu(3)-n3a
      do n4c=0,NumGlu(4)-n4a

         n1b = NumGlu(1)-n1a-n1c
         n3b = NumGlu(3)-n3a-n3c
         n4b = NumGlu(4)-n4a-n4c
         if( n1c.gt.0 .and. (n3c+n4c).gt.0  ) cycle
         if( n3c.gt.0 .and. (n1c+n4c).gt.0  ) cycle
         if( n4c.gt.0 .and. (n1c+n3c).gt.0  ) cycle

         rIn =n1a+n1c+1
         rOut=NumGlu(1)+NumGlu(2)+n3a
         Eps2 = cur_g_2f(Gluons(rIn:rOut),Quarks(2:3),(/1+n1b+NumGlu(2)+n3a,n1b,NumGlu(2),n3a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
         Eps2 = Eps2*PropFac1

            rIn =NumGlu(1)+NumGlu(2)+n3a+n3c+1
            rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
            Sca1 = cur_s_2s(Gluons(rIn:rOut),Scalar(4:4),(/n3b+n4a,n3b,n4a/) )! here was a bug
            PMom2(:) = Scalar(4)%Mom + SumMom(Gluons,rIn,rOut)! here was a bug
            if(n3b.ge.1 .or. n4a.ge.1) then
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Scalar(4)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Scalar(4)%Mass2).lt.PropCut ) cycle
               Sca1 = Sca1 * PropFac2
            endif

            if( n1c.gt.0 ) then
                rIn =n1a+1
                rOut=n1a+n1c
                EpsX(:) = cur_g(Gluons(rIn:rOut),1+n1c)
                if(n1c.gt.1) then
                    PMom4 = SumMom(Gluons,rIn,rOut)
                    PropFac4 = (0d0,-1d0)/sc_(PMom4,PMom4)
                    if( abs(sc_(PMom4,PMom4)).lt.PropCut ) cycle
                    EpsX(:) = EpsX(:)*PropFac4
                endif
                Sca0 = cggs(EpsX,Eps2) * Sca1
            elseif( n3c.gt.0) then
                rIn =NumGlu(1)+NumGlu(2)+n3a+1
                rOut=NumGlu(1)+NumGlu(2)+n3a+n3c
                EpsX(:) = cur_g(Gluons(rIn:rOut),1+n3c)
                if(n3c.gt.1) then
                    PMom4 = SumMom(Gluons,rIn,rOut)
                    PropFac4 = (0d0,-1d0)/sc_(PMom4,PMom4)
                    if( abs(sc_(PMom4,PMom4)).lt.PropCut ) cycle
                    EpsX(:) = EpsX(:)*PropFac4
                endif
                Sca0 = cggs(EpsX,Eps2) * Sca1
            elseif( n4c.gt.0 ) then
                rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+n4c
                EpsX(:) = cur_g(Gluons(rIn:rOut),1+n4c)
                if(n4c.gt.1) then
                    PMom4 = SumMom(Gluons,rIn,rOut)
                    PropFac4 = (0d0,-1d0)/sc_(PMom4,PMom4)
                    if( abs(sc_(PMom4,PMom4)).lt.PropCut ) cycle
                    EpsX(:) = EpsX(:)*PropFac4
                endif
                Sca0 = cgsg(Eps2,EpsX) * Sca1
            else
                if( Scalar(4)%PartType.gt.0 ) then
                    Sca0 = cgs(Eps2,PMom1,PMom1+PMom2) * Sca1
                else
                    Sca0 = cgbs(Eps2,PMom1,PMom1+PMom2) * Sca1
                endif
            endif

            PMom1 = Quarks(2)%Mom+Quarks(3)%Mom+Scalar(4)%Mom+SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+n4c)
            if(n1a.ge.1 .or. n4b.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Scalar(4)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Scalar(4)%Mass2).lt.PropCut ) cycle
               Sca0 = Sca0 * PropFac1
            endif


            TmpScalar(1)%Mom  => PMom1(:)
            TmpScalar(1)%Pol  => Sca0
            TmpScalar(1)%Mass => Scalar(4)%Mass
            TmpScalar(1)%Mass2=> Scalar(4)%Mass2
            TmpExtRef = -1
            TmpScalar(1)%ExtRef => TmpExtRef
            TmpScalar(1)%PartType => Scalar(4)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+n4c+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp = cur_s_2s(TmpGluons(1:counter-1),TmpScalar(1:1),(/counter-1,n1a,n4b/) )
            Res = Res + tmp

! print *, "bg res",tmp
! EpsX(1:Dv) = vbqq(dv, Quarks(3)%Pol,Quarks(2)%Pol )
! print *, "ms res",Scalar(4)%Pol(1) /dsqrt(2d0)/((Quarks(2)%Mom(1:Dv)+Quarks(3)%Mom(1:Dv)).ndot.(Quarks(2)%Mom(1:Dv)+Quarks(3)%Mom(1:Dv))) & 
!                  * (  (EpsX).ndot.( Scalar(4)%Mom+PMom1 )   )

      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

return
END FUNCTION





FUNCTION cur_f_ffss(Gluons,Scalars,Quarks,NumGlu) result(res)           ! Quarks(:) does not include the OFF-shell quark
implicit none                                                           ! checked gauge inv. for one gluon
integer :: NumGlu(0:4)
type(PtrToParticle) :: Gluons(1:),Quarks(2:2),Scalars(3:4)
integer,target :: TmpExtRef
complex(8) :: res(1:Ds),tmp(1:Ds)
complex(8) :: ubar1(1:Ds)
complex(8),target :: ubar0(1:Ds)
complex(8) :: eps1(1:Dv)
complex(8) :: eps2(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(4)),TmpQuark(1:1)
complex(8) :: PropFac1,PropFac2
complex(8),target :: pmom1(1:Dv)
complex(8) :: pmom2(1:Dv)
integer :: n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b
integer :: rIn,rOut,i,counter


!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4).ne.0 ) print *, "wrong number of gluons in cur_f_ffss"
!DEC$ ENDIF

   Res(:)=(0d0,0d0)
      do n2a=0,NumGlu(2)
      do n4a=0,NumGlu(4)
         n2b = NumGlu(2)-n2a
         n4b = NumGlu(4)-n4a

         rIn =NumGlu(1)+n2a+1
         rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
         Eps2 = cur_g_2s(Gluons(rIn:rOut),Scalars(3:4),(/1+n2b+NumGlu(3)+n4a,n2b,NumGlu(3),n4a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Scalars(3)%Mom + Scalars(4)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
         Eps2 = Eps2*PropFac1
         do n1a=0,NumGlu(1)
            n1b = NumGlu(1)-n1a
            rIn =n1a+1
            rOut=NumGlu(1)+n2a
            ubar1(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(2:2),-Quarks(2)%PartType,(/n2a+n1b,n1b,n2a/) )
            if(n1b.ge.1 .or. n2a.ge.1) then
               PMom2(:) = Quarks(2)%Mom + SumMom(Gluons,rIn,rOut)  ! can be moved outside the n1a-loop
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(2)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Quarks(2)%Mass2).lt.PropCut ) cycle
               if( Quarks(2)%PartType.lt.0 ) then
                  ubar1(:) = (-spi2_(PMom2,ubar1)+Quarks(2)%Mass*ubar1(:))*PropFac2
               else
                  ubar1(:) = (+spb2_(ubar1,PMom2)+Quarks(2)%Mass*ubar1(:))*PropFac2
               endif
            endif
            if( Quarks(2)%PartType.lt.0 ) then
               ubar0(:) = vbqg(ubar1,eps2)       ! re-checked
            else
               ubar0(:) = vqg(ubar1,eps2)        ! re-checked
            endif

            PMom1 = Quarks(2)%Mom+Scalars(3)%Mom+Scalars(4)%Mom+SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a)
            if(n1a.ge.1 .or. n4b.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(2)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Quarks(2)%Mass2).lt.PropCut ) cycle
               if( Quarks(2)%PartType.lt.0 ) then
                  ubar0(:) = (-spi2_(PMom1,ubar0)+Quarks(2)%Mass*ubar0(:))*PropFac1
               else
                  ubar0(:) = (+spb2_(ubar0,PMom1)+Quarks(2)%Mass*ubar0(:))*PropFac1
               endif
            endif

            TmpQuark(1)%Mom  => PMom1(:)
            TmpQuark(1)%Pol  => ubar0(:)
            TmpQuark(1)%Mass => Quarks(2)%Mass
            TmpQuark(1)%Mass2=> Quarks(2)%Mass2
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpQuark(1)%PartType => Quarks(2)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),-TmpQuark(1)%PartType,(/counter-1,n1a,n4b/) )
            Res(:) = Res(:) + tmp(:)
         enddo
      enddo
      enddo


return
END FUNCTION






FUNCTION cur_f_fssf(Gluons,Scalars,Quarks,NumGlu) result(res)           ! Quarks(:) does not include the OFF-shell quark
implicit none                                                           ! checked gauge inv. for one gluon
integer :: NumGlu(0:4)
type(PtrToParticle) :: Gluons(1:),Quarks(4:4),Scalars(2:3)
integer,target :: TmpExtRef
complex(8) :: res(1:Ds),tmp(1:Ds)
complex(8) :: ubar1(1:Ds)
complex(8),target :: ubar0(1:Ds)
complex(8) :: eps1(1:Dv)
complex(8) :: eps2(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(4)),TmpQuark(1:1)
complex(8) :: PropFac1,PropFac2
complex(8),target :: pmom1(1:Dv)
complex(8) :: pmom2(1:Dv)
integer :: n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b
integer :: rIn,rOut,i,counter


!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4).ne.0 ) print *, "wrong number of gluons in cur_f_fssf"
!DEC$ ENDIF

      Res(:) = (0d0,0d0)
      do n1a=0,NumGlu(1)
      do n3b=0,NumGlu(3)
         n1b = NumGlu(1)-n1a
         n3a = NumGlu(3)-n3b
         rIn =n1a+1
         rOut=NumGlu(1)+NumGlu(2)+n3a
         Eps2 = cur_g_2s(Gluons(rIn:rOut),Scalars(2:3),(/1+n1b+NumGlu(2)+n3a,n1b,NumGlu(2),n3a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Scalars(2)%Mom + Scalars(3)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
         Eps2 = Eps2*PropFac1

         do n4a=0,NumGlu(4)
            n4b = NumGlu(4)-n4a
            rIn =NumGlu(1)+NumGlu(2)+n3a+1
            rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
            ubar1(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(4:4),-Quarks(4)%PartType,(/NumGlu(3)+n4a-n3a,n3b,n4a/) )
            if(n3b.ge.1 .or. n4a.ge.1) then
               PMom2(:) = Quarks(4)%Mom + SumMom(Gluons,rIn,rOut)
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(4)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Quarks(4)%Mass2).lt.PropCut ) cycle
               if( Quarks(4)%PartType.lt.0 ) then
                  ubar1(:) = (-spi2_(PMom2,ubar1)+Quarks(4)%Mass*ubar1(:))*PropFac2
               else
                  ubar1(:) = (+spb2_(ubar1,PMom2)+Quarks(4)%Mass*ubar1(:))*PropFac2
               endif
            endif
            if( Quarks(4)%PartType.lt.0 ) then
               ubar0(:) = vgbq(eps2,ubar1)   !! changed from vqg(ubar1,eps2)       ! re-checked
            else
               ubar0(:) = vgq(eps2,ubar1)   !! changed from vbqg(ubar1,eps2)       ! re-checked
            endif

            PMom1 = Scalars(2)%Mom+Scalars(3)%Mom+Quarks(4)%Mom+SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a)
            if(n1a.ge.1 .or. n4b.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(4)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Quarks(4)%Mass2).lt.PropCut ) cycle
               if( Quarks(4)%PartType.lt.0 ) then
                  ubar0(:) = (-spi2_(PMom1,ubar0)+Quarks(4)%Mass*ubar0(:))*PropFac1
               else
                  ubar0(:) = (+spb2_(ubar0,PMom1)+Quarks(4)%Mass*ubar0(:))*PropFac1
               endif
            endif
            TmpQuark(1)%Mom  => PMom1(:)
            TmpQuark(1)%Pol  => ubar0(:)
            TmpQuark(1)%Mass => Quarks(4)%Mass
            TmpQuark(1)%Mass2=> Quarks(4)%Mass2
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpQuark(1)%PartType => Quarks(4)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),-TmpQuark(1)%PartType,(/counter-1,n1a,n4b/) )
            Res(:) = Res(:) + tmp(:)
         enddo
      enddo
      enddo


return
END FUNCTION






FUNCTION cur_f_fffssf(Gluons,Scalars,Quarks,NumGlu) result(res)           ! Quarks(:) does not include the OFF-shell quark
implicit none
integer :: NumGlu(0:6)
type(PtrToParticle) :: Gluons(1:),Quarks(2:4),Scalars(1:2)
integer,target :: TmpExtRef
complex(8) :: res(1:Ds),tmp(1:Ds)
complex(8) :: ubar1(1:Ds),u1(1:Ds)
complex(8),target :: ubar0(1:Ds)
complex(8) :: eps1(1:Dv)
complex(8) :: eps2(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(4)),TmpQuark(1:1)
complex(8) :: PropFac1,PropFac2
complex(8),target :: pmom1(1:Dv)
complex(8),target :: PMom2(1:Dv)
integer :: n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,n5a,n5b,n6a,n6b
integer :: rIn,rOut,i,counter


!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4)-NumGlu(5)-NumGlu(6).ne.0 ) print *, "wrong number of gluons in cur_f_fffssf"
!DEC$ ENDIF
    if( .not.(Quarks(2)%PartType.eq.-Quarks(3)%PartType .AND. abs(Quarks(2)%PartType).ne.abs(Quarks(4)%PartType)) ) call Error("this flavor combination is not yet implemented in cur_f_ffssf")


    Res(:) = (0d0,0d0)
      do n1a=0,NumGlu(1)
      do n5b=0,NumGlu(5)
         n1b = NumGlu(1)-n1a
         n5a = NumGlu(5)-n5b

         rIn =n1a+1
         rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a
         Eps2 = cur_g_ffss(Gluons(rIn:rOut),Scalars(1:2),Quarks(2:3),(/1+n1b+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a,n1b,NumGlu(2),NumGlu(3),NumGlu(4),n5a/)) 
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom + Scalars(1)%Mom + Scalars(2)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
         Eps2 = Eps2*PropFac1

         do n6a=0,NumGlu(6)
            n6b = NumGlu(6)-n6a
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a+1
            rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            ubar1(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(4:4),-Quarks(4)%PartType,(/n5b+n6a,n5b,n6a/) )
            if(n5b.ge.1 .or. n6a.ge.1) then
               PMom2(:) = SumMom(Gluons,rIn,rOut) + Quarks(4)%Mom
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(4)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Quarks(4)%Mass2).lt.PropCut ) cycle
               if( Quarks(4)%PartType.lt.0 ) then
                  ubar1(:) = (-spi2_(PMom2,ubar1)+Quarks(4)%Mass*ubar1(:))*PropFac2
               else
                  ubar1(:) = (+spb2_(ubar1,PMom2)+Quarks(4)%Mass*ubar1(:))*PropFac2
               endif
            endif
            if( Quarks(4)%PartType.lt.0 ) then
               ubar0(:) = vgbq(Eps2,ubar1)
            else
               ubar0(:) = vgq(Eps2,ubar1)
            endif

            rIn = n1a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            PMom1 = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom + Scalars(1)%Mom + Scalars(2)%Mom + Quarks(4)%Mom
            if(n1a.ge.1 .or. n6b.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(4)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Quarks(4)%Mass2).lt.PropCut ) cycle
               if( Quarks(4)%PartType.lt.0 ) then
                  ubar0(:) = (-spi2_(PMom1,ubar0)+Quarks(4)%Mass*ubar0(:))*PropFac1
               else
                  ubar0(:) = (+spb2_(ubar0,PMom1)+Quarks(4)%Mass*ubar0(:))*PropFac1
               endif
            endif
            TmpQuark(1)%Mom => PMom1(:)
            TmpQuark(1)%Pol => ubar0(:)
            TmpQuark(1)%Mass => Quarks(4)%Mass
            TmpQuark(1)%Mass2=> Quarks(4)%Mass2
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpQuark(1)%PartType => Quarks(4)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),-Quarks(4)%PartType,(/counter-1,n1a,n6b/) )

            Res(:) = Res(:) + tmp(:)
!             Res2(:) = Res2(:) + tmp(:)
!             print *, "2",tmp(:)
         enddo
      enddo
      enddo




      do n1a=0,NumGlu(1)
      do n3a=0,NumGlu(3)
      do n6a=0,NumGlu(6)
         n1b = NumGlu(1)-n1a
         n3b = NumGlu(3)-n3a
         n6b = NumGlu(6)-n6a

            counter=1
            rIn = n1a+1
            rOut= NumGlu(1)+NumGlu(2)+n3a
            Eps2 = cur_g_2f(Gluons(rIn:rOut),Quarks(2:3),(/1+n1b+NumGlu(2)+n3a,n1b,NumGlu(2),n3a/))
            PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
            PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
            if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
            Eps2 = Eps2*PropFac1

            rIn = NumGlu(1)+NumGlu(2)+n3a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            u1 = cur_f_fssf(Gluons(rIn:rOut),Scalars(1:2),Quarks(4:4),(/n3b+NumGlu(4)+NumGlu(5)+n6a,n3b,NumGlu(4),NumGlu(5),n6a/))
            PMom2  = SumMom(Gluons,rIn,rOut)  + Scalars(1)%Mom + Scalars(2)%Mom + Quarks(4)%Mom

                if( Quarks(4)%PartType.lt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = (-spi2_(PMom2,u1) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vgbq(eps2,u1)
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = (-spi2_(PMom2,ubar0) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                else 
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = ( spb2_(u1,PMom2) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vgq(eps2,u1)
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = ( spb2_(ubar0,PMom2) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                endif
                TmpQuark(1)%Mom  => PMom2(:)
                TmpQuark(1)%Pol  => ubar0(:)
                TmpQuark(1)%Mass => Quarks(4)%Mass
                TmpQuark(1)%Mass2=> Quarks(4)%Mass2

            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpQuark(1)%PartType => Quarks(4)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),-Quarks(4)%PartType,(/counter-1,n1a,n6b/))

            Res(:)  = Res(:) + Tmp(:)
!             Res4(:) = Res4(:) + Tmp(:)
!             print *, "4",tmp(:)
      enddo
      enddo
      enddo


return
END FUNCTION





FUNCTION cur_f_fssfff(Gluons,Scalars,Quarks,NumGlu) result(res)           ! Quarks(:) does not include the OFF-shell quark
implicit none
integer :: NumGlu(0:6)
type(PtrToParticle) :: Gluons(1:),Quarks(2:4),Scalars(1:2)
integer,target :: TmpExtRef
complex(8) :: res(1:Ds),tmp(1:Ds)
complex(8) :: ubar1(1:Ds),u1(1:Ds)
complex(8),target :: ubar0(1:Ds)
complex(8) :: eps1(1:Dv)
complex(8) :: eps2(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(4)),TmpQuark(1:1)
complex(8) :: PropFac1,PropFac2
complex(8),target :: pmom1(1:Dv)
complex(8),target :: PMom2(1:Dv)
integer :: n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,n5a,n5b,n6a,n6b
integer :: rIn,rOut,i,counter


!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4)-NumGlu(5)-NumGlu(6).ne.0 ) print *, "wrong number of gluons in cur_f_ssfff"
!DEC$ ENDIF
    if( .not.(Quarks(2)%PartType.eq.-Quarks(3)%PartType .AND. abs(Quarks(2)%PartType).ne.abs(Quarks(4)%PartType)) ) call Error("this flavor combination is not yet implemented in cur_f_ssfff")


    Res(:) = (0d0,0d0)
      do n1a=0,NumGlu(1)
      do n5b=0,NumGlu(5)
         n1b = NumGlu(1)-n1a
         n5a = NumGlu(5)-n5b

         rIn =n1a+1
         rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a
         Eps2 =  cur_g_ssff(Gluons(rIn:rOut),Scalars(1:2),Quarks(2:3),(/1+n1b+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a,n1b,NumGlu(2),NumGlu(3),NumGlu(4),n5a/)) 
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Scalars(1)%Mom + Scalars(2)%Mom + Quarks(2)%Mom + Quarks(3)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
         Eps2 = Eps2*PropFac1

         do n6a=0,NumGlu(6)
            n6b = NumGlu(6)-n6a
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a+1
            rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            ubar1(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(4:4),-Quarks(4)%PartType,(/n5b+n6a,n5b,n6a/) )
            if(n5b.ge.1 .or. n6a.ge.1) then
               PMom2(:) = SumMom(Gluons,rIn,rOut) + Quarks(4)%Mom
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(4)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Quarks(4)%Mass2).lt.PropCut ) cycle
               if( Quarks(4)%PartType.lt.0 ) then
                  ubar1(:) = (-spi2_(PMom2,ubar1)+Quarks(4)%Mass*ubar1(:))*PropFac2
               else
                  ubar1(:) = (+spb2_(ubar1,PMom2)+Quarks(4)%Mass*ubar1(:))*PropFac2
               endif
            endif
            if( Quarks(4)%PartType.lt.0 ) then
               ubar0(:) = vgbq(Eps2,ubar1)
            else
               ubar0(:) = vgq(Eps2,ubar1)
            endif

            rIn = n1a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            PMom1 = SumMom(Gluons,rIn,rOut) + Scalars(1)%Mom + Scalars(2)%Mom + Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom
            if(n1a.ge.1 .or. n6b.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(4)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Quarks(4)%Mass2).lt.PropCut ) cycle
               if( Quarks(4)%PartType.lt.0 ) then
                  ubar0(:) = (-spi2_(PMom1,ubar0)+Quarks(4)%Mass*ubar0(:))*PropFac1
               else
                  ubar0(:) = (+spb2_(ubar0,PMom1)+Quarks(4)%Mass*ubar0(:))*PropFac1
               endif
            endif
            TmpQuark(1)%Mom => PMom1(:)
            TmpQuark(1)%Pol => ubar0(:)
            TmpQuark(1)%Mass => Quarks(4)%Mass
            TmpQuark(1)%Mass2=> Quarks(4)%Mass2
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpQuark(1)%PartType => Quarks(4)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),-Quarks(4)%PartType,(/counter-1,n1a,n6b/) )

            Res(:) = Res(:) + tmp(:)
!             Res2(:) = Res2(:) + tmp(:)
!             print *, "2",tmp(:)
         enddo
      enddo
      enddo



      do n1a=0,NumGlu(1)
      do n3a=0,NumGlu(3)
      do n6a=0,NumGlu(6)
         n1b = NumGlu(1)-n1a
         n3b = NumGlu(3)-n3a
         n6b = NumGlu(6)-n6a

            counter=1
            rIn = n1a+1
            rOut= NumGlu(1)+NumGlu(2)+n3a
            Eps2 = cur_g_2s(Gluons(rIn:rOut),Scalars(1:2),(/1+n1b+NumGlu(2)+n3a,n1b,NumGlu(2),n3a/))
            PMom1(:) = SumMom(Gluons,rIn,rOut) + Scalars(1)%Mom + Scalars(2)%Mom
            PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
            if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
            Eps2 = Eps2*PropFac1

            rIn = NumGlu(1)+NumGlu(2)+n3a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            u1 = cur_f_4f(Gluons(rIn:rOut),Quarks(2:4),-Quarks(4)%PartType,(/n3b+NumGlu(4)+NumGlu(5)+n6a,n3b,NumGlu(4),NumGlu(5),n6a/),2)
            PMom2  = SumMom(Gluons,rIn,rOut)  + Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom

                if( Quarks(4)%PartType.lt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = (-spi2_(PMom2,u1) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vgbq(eps2,u1) 
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = (-spi2_(PMom2,ubar0) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                elseif( Quarks(4)%PartType.gt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = ( spb2_(u1,PMom2) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vgq(eps2,u1) 
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = ( spb2_(ubar0,PMom2) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                endif
                TmpQuark(1)%Mom  => PMom2(:)
                TmpQuark(1)%Pol  => ubar0(:)
                TmpQuark(1)%Mass => Quarks(4)%Mass
                TmpQuark(1)%Mass2=> Quarks(4)%Mass2
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpQuark(1)%PartType => Quarks(4)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),-Quarks(4)%PartType,(/counter-1,n1a,n6b/))

            Res(:)  = Res(:) + Tmp(:)
!             Res4(:) = Res4(:) + Tmp(:)
!             print *, "4",tmp(:)
      enddo
      enddo
      enddo


return
END FUNCTION









FUNCTION cur_g_2s(Gluons,Scalars,NumGlu) result(Res)           ! Gluons(:) does not include the OFF-shell gluon,  however NumGlu(0) is the number of all gluons
implicit none                                                  ! checked gauge invariance for up to 4 gluons
integer :: NumGlu(0:3),i,counter
type(PtrToParticle) :: Gluons(1:),Scalars(1:2)
integer :: rIn,rOut,n1a,n1b,n1c,n2a,n2b,n2c,n3a,n3b,n3c
integer,target :: TmpExtRef
complex(8) :: Res(1:Dv)
complex(8) :: sc1,sc2
complex(8),target :: Eps1(1:Dv)
complex(8) :: Eps2(1:Dv)
complex(8) :: EpsX(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(3)+1)
complex(8) :: PMom1(1:Dv),PMom2(1:Dv),PMom4(1:Dv)
complex(8),target :: PMom3(1:Dv)
complex(8) :: PropFac1,PropFac2,PropFac3,PropFac4
integer :: PartKey,HelKey,CurrKey,Hel_Tmp


!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-1-NumGlu(1)-NumGlu(2)-NumGlu(3).ne.0 ) print *, "wrong number of gluons in cur_g_2s"
!DEC$ ENDIF

   res = (0d0,0d0)
   do n1a=0,NumGlu(1)
   do n3a=0,NumGlu(3)
   do n2a=0,NumGlu(2)
   do n1c=0,NumGlu(1)-n1a
   do n2c=0,NumGlu(2)-n2a
   do n3c=0,NumGlu(3)-n3a


      if( n1c.gt.0 .and. (n2c+n3c).gt.0  ) cycle
      if( n2c.gt.0 .and. (n1c+n3c).gt.0  ) cycle
      if( n3c.gt.0 .and. (n1c+n2c).gt.0  ) cycle

      n1b=NumGlu(1)-n1a-n1c
      n2b=NumGlu(2)-n2a-n2c
      n3b=NumGlu(3)-n3a-n3c

      ! Fer1
      rIn=n1a+n1c+1
      rOut=NumGlu(1)+n2a
      PMom1(:) = Scalars(1)%Mom(:)+ SumMom(Gluons,rIn,rOut)
      sc1 = cur_s_2s(Gluons(rIn:rOut),Scalars(1:1),(/n1b+n2a,n1b,n2a/))
      if(n1b.ge.1 .or. n2a.ge.1) then
         PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Scalars(1)%Mass2)
         if( abs(sc_(PMom1,PMom1)-Scalars(1)%Mass2).lt.PropCut ) cycle
         sc1 = sc1*PropFac1
      endif

      ! Fer2
      rIn=NumGlu(1)+n2a+n2c+1
      rOut=NumGlu(1)+NumGlu(2)+n3a
      PMom2(:) = Scalars(2)%Mom(:)+ SumMom(Gluons,rIn,rOut)
      sc2 = cur_s_2s(Gluons(rIn:rOut),Scalars(2:2),(/n2b+n3a,n2b,n3a/))
      if(n2b.ge.1 .or. n3a.ge.1) then
         PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Scalars(2)%Mass2)
         if( abs(sc_(PMom2,PMom2)-Scalars(2)%Mass2).lt.PropCut ) cycle
         sc2 = sc2*PropFac2
      endif

      if( n1c.gt.0 ) then
          rIn =n1a+1
          rOut=n1a+n1c
          EpsX(:) = cur_g(Gluons(rIn:rOut),1+n1c)
          if(n1c.gt.1) then
              PMom4 = SumMom(Gluons,rIn,rOut)
              PropFac4 = (0d0,-1d0)/sc_(PMom4,PMom4)
              if( abs(sc_(PMom4,PMom4)).lt.PropCut ) cycle
              EpsX(:) = EpsX(:)*PropFac4
          endif
          Eps1(:) = vggss(Dv,EpsX) * sc1*sc2 

      elseif( n2c.gt.0) then
          rIn =NumGlu(1)+n2a+1
          rOut=NumGlu(1)+n2a+n2c
          EpsX(:) = cur_g(Gluons(rIn:rOut),1+n2c)
          if(n2c.gt.1) then
              PMom4 = SumMom(Gluons,rIn,rOut)
              PropFac4 = (0d0,-1d0)/sc_(PMom4,PMom4)
              if( abs(sc_(PMom4,PMom4)).lt.PropCut ) cycle
              EpsX(:) = EpsX(:)*PropFac4
          endif
          Eps1(:) = vgsgs(Dv,EpsX) * sc1*sc2

      elseif( n3c.gt.0 ) then
          rIn =NumGlu(1)+NumGlu(2)+n3a+1
          rOut=NumGlu(1)+NumGlu(2)+n3a+n3c
          EpsX(:) = cur_g(Gluons(rIn:rOut),1+n3c)
          if(n3c.gt.1) then
              PMom4 = SumMom(Gluons,rIn,rOut)
              PropFac4 = (0d0,-1d0)/sc_(PMom4,PMom4)
              if( abs(sc_(PMom4,PMom4)).lt.PropCut ) cycle
              EpsX(:) = EpsX(:)*PropFac4
          endif
          Eps1(:) = vgssg(Dv,EpsX) * sc1*sc2
      else
          if( Scalars(1)%PartType.gt.0 ) then
             Eps1(:) = vbss(Dv,PMom1(:),PMom2(:))*sc1*sc2
          else
             Eps1(:) = vsbs(Dv,PMom1(:),PMom2(:))*sc1*sc2
          endif
      endif

      PMom3(:) = Scalars(1)%Mom(:)+Scalars(2)%Mom(:) + SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+n3a+n3c)! here was a bug

      counter=1
      rIn =1
      rOut=n1a
      do i=rIn,rOut
         call CopyParticlePtr(Gluons(i),TmpGluons(counter))
         counter=counter+1
      enddo
      TmpGluons(counter)%Mom => PMom3(:)
      TmpGluons(counter)%Pol => Eps1(:)
      TmpExtRef = -1
      TmpGluons(counter)%ExtRef => TmpExtRef
      counter=counter+1
      rIn =NumGlu(1)+NumGlu(2)+n3a+n3c+1
      rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)
      do i=rIn,rOut
         call CopyParticlePtr(Gluons(i),TmpGluons(counter))
         counter=counter+1
      enddo
      Eps2(:) = cur_g(TmpGluons(1:counter-1),1+n1a+n3b+1)

      if(n1a.ge.1 .or. n3b.ge.1) then
         PropFac3 = (0d0,-1d0)/sc_(PMom3,PMom3)
         if( abs(sc_(PMom3,PMom3)).lt.PropCut ) cycle
         Eps2(:) = Eps2(:)*PropFac3
      endif

      Res(:) = Res(:) + Eps2(:)
   enddo
   enddo
   enddo
   enddo
   enddo
   enddo


! if(  numglu(0).eq.4 ) then
! print *, "cur_g_2s",Res(1:4)
! pause
! endif


return
END FUNCTION














FUNCTION cur_g_ssff(Gluons,Scalars,Quarks,NumGlu) result(res)           ! Gluons(:) does not include the OFF-shell gluon,  however NumGlu is the number of all gluons
implicit none                                                           ! checked gauge inv. for one gluon                                                        
integer,intent(in) :: NumGlu(0:5)
type(PtrToParticle) :: Gluons(1:),Scalars(1:2),Quarks(3:4)
integer :: na,nb,nc,nd,ne,nf,ng,nh,ni,nj,nk
integer :: rIn,rOut
integer :: counter,i
complex(8) :: res(Dv)
type(PtrToParticle) :: TmpGluons(1:2+NumGlu(1)+NumGlu(3)+NumGlu(5))
complex(8),target :: TmpMom1(1:Dv),TmpMom2(1:Dv)
integer,target :: TmpExtRef1,TmpExtRef2
complex(8),target :: Eps1(1:Dv)
complex(8),target :: Eps2(1:Dv)
complex(8) :: Eps3(1:Dv)
complex(8) :: Sca1,Sca2
complex(8) :: u1(1:Ds),ubar2(1:Ds)
complex(8) :: PropFac1,PropFac2,PropFac3,PropFac4
complex(8) :: PMom1(1:Dv)
complex(8) :: PMom2(1:Dv)
complex(8) :: PMom3(1:Dv)
complex(8) :: PMom4(1:Dv)

!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-1-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4)-NumGlu(5).ne.0 ) print *, "wrong number of gluons in cur_g_ssff"
!DEC$ ENDIF

      res = (0d0,0d0)
!        type (1)
         do na=0,NumGlu(1)
         do nc=0,NumGlu(2)
         do ne=0,NumGlu(5)
            nb=NumGlu(1)-na
            nd=NumGlu(2)-nc
            nf=NumGlu(5)-ne

            rIn = NumGlu(1)+nc+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne
            Sca1 = cur_s_ssff(Gluons(rIn:rOut),Scalars(2:2),Quarks(3:4),(/nd+NumGlu(3)+NumGlu(4)+ne,nd,NumGlu(3),NumGlu(4),ne/))
            PMom2  = SumMom(Gluons,rIn,rOut)  + Scalars(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom
            PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Scalars(1)%Mass2)
            if( abs(sc_(PMom2,PMom2) - Scalars(1)%Mass2).lt.PropCut ) cycle
            Sca1 = Sca1*PropFac2

            rIn = na+1
            rOut= NumGlu(1)+nc
            Sca2 = cur_s_2s(Gluons(rIn:rOut),Scalars(1:1),(/nb+nc,nb,nc/))
            PMom3 = SumMom(Gluons,rIn,rOut) + Scalars(1)%Mom
            if( nb.ge.1 .or. nc.ge.1 ) then
               PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Scalars(1)%Mass2)
               if( abs(sc_(PMom3,PMom3) - Scalars(1)%Mass2).lt.PropCut ) cycle
               Sca2 = Sca2*PropFac3
            endif

            if( Scalars(1)%PartType.gt.0 ) then
                Eps1(:) = vbss(Dv,PMom3(:),PMom2(:)) * Sca1*Sca2! here was a bug: PMom2 <-->PMom3
            else
                Eps1(:) = vsbs(Dv,PMom3(:),PMom2(:)) * Sca1*Sca2! here was a bug: PMom2 <-->PMom3
            endif
            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpMom1(:) = PMom2(:)+PMom3(:)
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps1(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps2(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+1)

            if( na.ge.1 .or. nf.ge.1 ) then
               PropFac1 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
               if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
               Eps2 = Eps2*PropFac1
            endif

            Res = Res + Eps2
         enddo
         enddo
         enddo


!        type (2)
         do na=0,NumGlu(1)
         do nc=0,NumGlu(4)
         do ne=0,NumGlu(5)
            nb=NumGlu(1)-na
            nd=NumGlu(4)-nc
            nf=NumGlu(5)-ne

            rIn = na+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+nc
            u1 = cur_f_fssf(Gluons(rIn:rOut),Scalars(1:2),Quarks(3:3),(/nb+NumGlu(2)+NumGlu(3)+nc,nb,NumGlu(2),NumGlu(3),nc/))
            PMom2 =  SumMom(Gluons,rIn,rOut) + Scalars(1)%Mom + Scalars(2)%Mom + Quarks(3)%Mom
            PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
            if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
            if( Quarks(4)%PartType.lt.0 ) then
              u1 = (+spb2_(u1,PMom2) + Quarks(4)%Mass*u1 )*PropFac2
            else
              u1 = (-spi2_(PMom2,u1) + Quarks(4)%Mass*u1 )*PropFac2
            endif
            rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+nc+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne
            ubar2 = cur_f_2f(Gluons(rIn:rOut),Quarks(4:4),-Quarks(4)%PartType,(/nd+ne,nd,ne/))
            PMom3  = SumMom(Gluons,rIn,rOut) + Quarks(4)%Mom
            if( nd.ge.1 .or. ne.ge.1 ) then
               PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Quarks(4)%Mass2)
               if( abs(sc_(PMom3,PMom3) - Quarks(4)%Mass2).lt.PropCut ) cycle
               if( Quarks(4)%PartType.lt.0 ) then
                  ubar2 = (-spi2_(PMom3,ubar2) + Quarks(4)%Mass*ubar2 )*PropFac3
               else
                  ubar2 = (+spb2_(ubar2,PMom3) + Quarks(4)%Mass*ubar2 )*PropFac3
               endif
            endif

            if( Quarks(4)%PartType.lt.0 ) then
              Eps1 = +vbqq(Dv,u1,ubar2)
            else
              Eps1 = -vbqq(Dv,ubar2,u1)
            endif


            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpMom1(:) = PMom2(:)+PMom3(:)
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps1(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps2(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+1)

            if( na.ge.1 .or. nf.ge.1 ) then
               PropFac1 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
               if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
               Eps2 = Eps2*PropFac1
            endif
            Res = Res + Eps2
         enddo
         enddo
         enddo



!        type(3)
         do na=0,NumGlu(1)
         do nc=0,NumGlu(2)
         do ne=0,NumGlu(3)
         do nf=0,NumGlu(3)-ne
         do nh=0,NumGlu(4)
         do nj=0,NumGlu(5)
            nb=NumGlu(1)-na
            nd=NumGlu(2)-nc
            ng=NumGlu(3)-ne-nf
            ni=NumGlu(4)-nh
            nk=NumGlu(5)-nj

            rIn = na+1
            rOut= NumGlu(1)+nc
            Sca2 = cur_s_2s(Gluons(rIn:rOut),Scalars(1:1),(/nb+nc,nb,nc/))
            PMom1 = SumMom(Gluons,rIn,rOut) + Scalars(1)%Mom
            if( nb.ge.1 .or. nc.ge.1 ) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1) - Scalars(1)%Mass2)
               if( abs(sc_(PMom1,PMom1) - Scalars(1)%Mass2).lt.PropCut ) cycle
               Sca2 = Sca2*PropFac1
            endif

            rIn = NumGlu(1)+nc+1
            rOut= NumGlu(1)+NumGlu(2)+ne
            Sca1 = cur_s_2s(Gluons(rIn:rOut),Scalars(2:2),(/nd+ne,nd,ne/))
            PMom2 = SumMom(Gluons,rIn,rOut) + Scalars(2)%Mom
            if( nd.ge.1 .or. ne.ge.1 ) then
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Scalars(2)%Mass2)
               if( abs(sc_(PMom2,PMom2) - Scalars(2)%Mass2).lt.PropCut ) cycle
               Sca1 = Sca1*PropFac2
            endif

            if( Scalars(1)%PartType.gt.0 ) then
                Eps1(:) = vbss(Dv,PMom1(:),PMom2(:)) * Sca1*Sca2! checked
            else
                Eps1(:) = vsbs(Dv,PMom1(:),PMom2(:)) * Sca1*Sca2! checked
            endif
            TmpMom1 = PMom1 + PMom2
            PropFac3 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
            if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
            Eps1 = Eps1*PropFac3


            rIn = NumGlu(1)+NumGlu(2)+ne+nf+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+nh
            u1 = cur_f_2f(Gluons(rIn:rOut),Quarks(3:3),-Quarks(3)%PartType,(/ng+nh,ng,nh/))
            PMom3 = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom
            if( ng.ge.1 .or. nh.ge.1 ) then
               PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Quarks(3)%Mass2)
               if( abs(sc_(PMom3,PMom3) - Quarks(3)%Mass2).lt.PropCut ) cycle
               if( Quarks(3)%PartType.lt.0 ) then
                  u1 = (-spi2_(PMom3,u1) + Quarks(3)%Mass*u1 )*PropFac3
               else
                  u1 = (+spb2_(u1,PMom3) + Quarks(3)%Mass*u1 )*PropFac3
               endif
            endif
            rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+nh+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+nj
            ubar2 = cur_f_2f(Gluons(rIn:rOut),Quarks(4:4),-Quarks(4)%PartType,(/ni+nj,ni,nj/))
            PMom4 = SumMom(Gluons,rIn,rOut) + Quarks(4)%Mom
            if( ni.ge.1 .or. nj.ge.1 ) then
               PropFac4 = (0d0,1d0)/(sc_(PMom4,PMom4) - Quarks(4)%Mass2)
               if( abs(sc_(PMom4,PMom4) - Quarks(4)%Mass2).lt.PropCut ) cycle
               if( Quarks(4)%PartType.lt.0 ) then
                  ubar2 = (-spi2_(PMom4,ubar2) + Quarks(4)%Mass*ubar2 )*PropFac4
               else
                  ubar2 = (+spb2_(ubar2,PMom4) + Quarks(4)%Mass*ubar2 )*PropFac4
               endif
            endif

            if( Quarks(4)%PartType.lt.0 ) then
               Eps2 = +vbqq(Dv,u1,ubar2)       ! re-checked
            else
               Eps2 = -vbqq(Dv,ubar2,u1)       ! re-checked
            endif
            TmpMom2 = PMom3 + PMom4
            PropFac1 = (0d0,-1d0)/sc_(TmpMom2,TmpMom2)
            if( abs(sc_(TmpMom2,TmpMom2)).lt.PropCut ) cycle
            Eps2 = Eps2*PropFac1


            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps1(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+ne+1
            rOut=NumGlu(1)+NumGlu(2)+ne+nf
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpGluons(counter)%Mom => TmpMom2(:)
            TmpGluons(counter)%Pol => Eps2(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+nj+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps3(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+nk+2)
            Res = Res + Eps3
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo



!        type (4)
         do na=0,NumGlu(1)
         do nc=0,NumGlu(2)
         do ne=0,NumGlu(3)
         do nj=0,NumGlu(5)
            nb=NumGlu(1)-na
            nd=NumGlu(2)-nc
            nf=NumGlu(3)-ne
            nk=NumGlu(5)-nj

            rIn = na+1
            rOut= NumGlu(1)+nc
            Sca2 = cur_s_2s(Gluons(rIn:rOut),Scalars(1:1),(/nb+nc,nb,nc/))
            PMom1 = SumMom(Gluons,rIn,rOut) + Scalars(1)%Mom
            if( nb.ge.1 .or. nc.ge.1 ) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1) - Scalars(1)%Mass2)
               if( abs(sc_(PMom1,PMom1) - Scalars(1)%Mass2).lt.PropCut ) cycle
               Sca2 = Sca2*PropFac1
            endif

            rIn = NumGlu(1)+nc+1
            rOut= NumGlu(1)+NumGlu(2)+ne
            Sca1 = cur_s_2s(Gluons(rIn:rOut),Scalars(2:2),(/nd+ne,nd,ne/))
            PMom2 = SumMom(Gluons,rIn,rOut) + Scalars(2)%Mom
            if( nd.ge.1 .or. ne.ge.1 ) then
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Scalars(2)%Mass2)
               if( abs(sc_(PMom2,PMom2) - Scalars(2)%Mass2).lt.PropCut ) cycle
               Sca1 = Sca1*PropFac2
            endif

            rIn = NumGlu(1)+NumGlu(2)+ne+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+nj
            Eps1(:) = cur_g_2f(Gluons(rIn:rOut),Quarks(3:4),(/1+nf+NumGlu(4)+nj,nf,NumGlu(4),nj/))
            PMom3 = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom + Quarks(4)%Mom
            PropFac3 = (0d0,-1d0)/sc_(PMom3,PMom3)! here was a bug: minus sign was missing
            if( abs(sc_(PMom3,PMom3)).lt.PropCut ) cycle
            Eps1 = Eps1*PropFac3


            Eps2(:) = vgssg(Dv,Eps1) * Sca1*Sca2


            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpMom1(:) = PMom1(:)+PMom2(:)+PMom3(:)
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps2(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+nj+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps3(:) = cur_g(TmpGluons(1:counter-1),1+na+nk+1)

            if( na.ge.1 .or. nk.ge.1 ) then
               PropFac1 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
               if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
               Eps3 = Eps3*PropFac1
            endif

            Res = Res + Eps3
         enddo
         enddo
         enddo
         enddo

return
END FUNCTION




FUNCTION cur_g_ffss(Gluons,Scalars,Quarks,NumGlu) result(res)           ! Gluons(:) does not include the OFF-shell gluon,  however NumGlu is the number of all gluons
implicit none                                                           ! checked gauge inv. for one gluon
integer,intent(in) :: NumGlu(0:5)
type(PtrToParticle) :: Gluons(1:),Quarks(1:2),Scalars(3:4)
integer :: na,nb,nc,nd,ne,nf,ng,nh,ni,nj,nk
integer :: rIn,rOut
integer :: counter,i
complex(8) :: res(Dv)
type(PtrToParticle) :: TmpGluons(1:2+NumGlu(1)+NumGlu(3)+NumGlu(5))
complex(8),target :: TmpMom1(1:Dv),TmpMom2(1:Dv)
integer,target :: TmpExtRef1,TmpExtRef2
complex(8),target :: Eps1(1:Dv)
complex(8),target :: Eps2(1:Dv)
complex(8) :: Eps3(1:Dv)
complex(8) :: u1(1:Ds),ubar2(1:Ds)
complex(8) :: Sca1,Sca2
complex(8) :: PropFac1,PropFac2,PropFac3,PropFac4
complex(8) :: PMom1(1:Dv)
complex(8) :: PMom2(1:Dv)
complex(8) :: PMom3(1:Dv)
complex(8) :: PMom4(1:Dv)

!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-1-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4)-NumGlu(5).ne.0 ) print *, "wrong number of gluons in cur_g_ffss"
!DEC$ ENDIF

      res = (0d0,0d0)
!        type (1)
         do na=0,NumGlu(1)
         do nc=0,NumGlu(2)
         do ne=0,NumGlu(5)
            nb=NumGlu(1)-na
            nd=NumGlu(2)-nc
            nf=NumGlu(5)-ne

            rIn = NumGlu(1)+nc+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne
            u1 = cur_f_ffss(Gluons(rIn:rOut),Scalars(3:4),Quarks(2:2),(/nd+NumGlu(3)+NumGlu(4)+ne,nd,NumGlu(3),NumGlu(4),ne/))
            PMom2  = SumMom(Gluons,rIn,rOut)  + Quarks(2)%Mom + Scalars(3)%Mom + Scalars(4)%Mom
            PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(1)%Mass2)
            if( abs(sc_(PMom2,PMom2) - Quarks(1)%Mass2).lt.PropCut ) cycle
            if( Quarks(1)%PartType.lt.0 ) then
               u1 = ( spb2_(u1,PMom2) + Quarks(1)%Mass*u1 )*PropFac2
            else
               u1 = (-spi2_(PMom2,u1) + Quarks(1)%Mass*u1 )*PropFac2
            endif

            rIn = na+1
            rOut= NumGlu(1)+nc
            ubar2 = cur_f_2f(Gluons(rIn:rOut),Quarks(1:1),-Quarks(1)%PartType,(/nb+nc,nb,nc/))
            PMom3 = SumMom(Gluons,rIn,rOut) + Quarks(1)%Mom
            if( nb.ge.1 .or. nc.ge.1 ) then
               PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Quarks(1)%Mass2)
               if( abs(sc_(PMom3,PMom3) - Quarks(1)%Mass2).lt.PropCut ) cycle
               if( Quarks(1)%PartType.lt.0 ) then
                  ubar2 = (-spi2_(PMom3,ubar2) + Quarks(1)%Mass*ubar2 )*PropFac3
               else
                  ubar2 = (+spb2_(ubar2,PMom3) + Quarks(1)%Mass*ubar2 )*PropFac3
               endif
            endif

            if( Quarks(1)%PartType.lt.0 ) then
                Eps1 = -vbqq(Dv,u1,ubar2)       ! re-checked
            else
                Eps1 = +vbqq(Dv,ubar2,u1)       ! re-checked
            endif

            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpMom1(:) = PMom2(:)+PMom3(:)
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps1(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps2(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+1)

            if( na.ge.1 .or. nf.ge.1 ) then
               PropFac1 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
               if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
               Eps2 = Eps2*PropFac1
            endif

            Res = Res + Eps2
         enddo
         enddo
         enddo



!        type (2)
         do na=0,NumGlu(1)
         do nc=0,NumGlu(4)
         do ne=0,NumGlu(5)
            nb=NumGlu(1)-na
            nd=NumGlu(4)-nc
            nf=NumGlu(5)-ne

            rIn = na+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+nc
            Sca1 = cur_s_sffs(Gluons(rIn:rOut),Scalars(3:3),Quarks(1:2),(/nb+NumGlu(2)+NumGlu(3)+nc,nb,NumGlu(2),NumGlu(3),nc/))
            PMom2 =  SumMom(Gluons,rIn,rOut) + Quarks(1)%Mom + Quarks(2)%Mom + Scalars(3)%Mom
            PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Scalars(3)%Mass2)
            if( abs(sc_(PMom2,PMom2) - Scalars(4)%Mass2).lt.PropCut ) cycle
            Sca1 = Sca1 * PropFac2

            rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+nc+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne
            Sca2 = cur_s_2s(Gluons(rIn:rOut),Scalars(4:4),(/nd+ne,nd,ne/))
            PMom3  = SumMom(Gluons,rIn,rOut) + Scalars(4)%Mom
            if( nd.ge.1 .or. ne.ge.1 ) then
               PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Scalars(4)%Mass2)
               if( abs(sc_(PMom3,PMom3) - Scalars(4)%Mass2).lt.PropCut ) cycle
               Sca2 = Sca2 * PropFac3
            endif

            if( Scalars(3)%PartType.gt.0 ) then
                Eps1(:) = vbss(Dv,PMom2(:),PMom3(:)) * Sca1*Sca2! checked
            else
                Eps1(:) = vsbs(Dv,PMom2(:),PMom3(:)) * Sca1*Sca2! checked
            endif

            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpMom1(:) = PMom2(:)+PMom3(:)
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps1(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps2(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+1)

            if( na.ge.1 .or. nf.ge.1 ) then
               PropFac1 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
               if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
               Eps2 = Eps2*PropFac1
            endif

            Res = Res + Eps2
         enddo
         enddo
         enddo


!        type(3)
         do na=0,NumGlu(1)
         do nc=0,NumGlu(2)
         do ne=0,NumGlu(3)
         do nf=0,NumGlu(3)-ne
         do nh=0,NumGlu(4)
         do nj=0,NumGlu(5)

            nb=NumGlu(1)-na
            nd=NumGlu(2)-nc
            ng=NumGlu(3)-ne-nf
            ni=NumGlu(4)-nh
            nk=NumGlu(5)-nj

            rIn = na+1
            rOut= NumGlu(1)+nc
            ubar2 = cur_f_2f(Gluons(rIn:rOut),Quarks(1:1),-Quarks(1)%PartType,(/nb+nc,nb,nc/))
            PMom1 = SumMom(Gluons,rIn,rOut) + Quarks(1)%Mom
            if( nb.ge.1 .or. nc.ge.1 ) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1) - Quarks(1)%Mass2)
               if( abs(sc_(PMom1,PMom1) - Quarks(1)%Mass2).lt.PropCut ) cycle
               if( Quarks(1)%PartType.lt.0 ) then
                  ubar2 = (-spi2_(PMom1,ubar2) + Quarks(1)%Mass*ubar2 )*PropFac1
               else
                  ubar2 = (+spb2_(ubar2,PMom1) + Quarks(1)%Mass*ubar2 )*PropFac1
               endif
            endif

            rIn = NumGlu(1)+nc+1
            rOut= NumGlu(1)+NumGlu(2)+ne
            u1 = cur_f_2f(Gluons(rIn:rOut),Quarks(2:2),-Quarks(2)%PartType,(/nd+ne,nd,ne/))
            PMom2 = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom
            if( nd.ge.1 .or. ne.ge.1 ) then
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
               if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
               if( Quarks(2)%PartType.lt.0 ) then
                  u1 = (-spi2_(PMom2,u1) + Quarks(2)%Mass*u1 )*PropFac2
               else
                  u1 = (+spb2_(u1,PMom2) + Quarks(2)%Mass*u1 )*PropFac2
               endif
            endif

            if( Quarks(2)%PartType.lt.0 ) then
              Eps1 = +vbqq(Dv,ubar2,u1)       ! re-checked
            else
              Eps1 = -vbqq(Dv,u1,ubar2)       ! re-checked
            endif
            TmpMom1 = PMom1 + PMom2
            PropFac3 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
            if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
            Eps1 = Eps1*PropFac3



            rIn = NumGlu(1)+NumGlu(2)+ne+nf+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+nh
            Sca1 = cur_s_2s(Gluons(rIn:rOut),Scalars(3:3),(/ng+nh,ng,nh/))
            PMom3 = SumMom(Gluons,rIn,rOut) + Scalars(3)%Mom
            if( ng.ge.1 .or. nh.ge.1 ) then
               PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Scalars(3)%Mass2)
               if( abs(sc_(PMom3,PMom3) - Scalars(3)%Mass2).lt.PropCut ) cycle
               Sca1 = Sca1 * PropFac3
            endif

            rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+nh+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+nj
            Sca2 = cur_s_2s(Gluons(rIn:rOut),Scalars(4:4),(/ni+nj,ni,nj/))
            PMom4 = SumMom(Gluons,rIn,rOut) + Scalars(4)%Mom
            if( ni.ge.1 .or. nj.ge.1 ) then
               PropFac4 = (0d0,1d0)/(sc_(PMom4,PMom4) - Scalars(4)%Mass2)
               if( abs(sc_(PMom4,PMom4) - Scalars(4)%Mass2).lt.PropCut ) cycle
               Sca2 = Sca2 * PropFac4
            endif

            if( Scalars(3)%PartType.gt.0 ) then
                Eps2(:) = vbss(Dv,PMom3(:),PMom4(:)) * Sca1*Sca2
            else
                Eps2(:) = vsbs(Dv,PMom3(:),PMom4(:)) * Sca1*Sca2
            endif
            TmpMom2 = PMom3 + PMom4
            PropFac3 = (0d0,-1d0)/sc_(TmpMom2,TmpMom2)
            if( abs(sc_(TmpMom2,TmpMom2)).lt.PropCut ) cycle
            Eps2 = Eps2*PropFac3


            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps1(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+ne+1
            rOut=NumGlu(1)+NumGlu(2)+ne+nf
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpGluons(counter)%Mom => TmpMom2(:)
            TmpGluons(counter)%Pol => Eps2(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+nj+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps3(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+nk+2)
            Res = Res + Eps3
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo



!        type(4)
         do na=0,NumGlu(1)
         do ne=0,NumGlu(3)
         do nh=0,NumGlu(4)
         do nj=0,NumGlu(5)

            nb=NumGlu(1)-na
            nf=NumGlu(3)-ne
            ni=NumGlu(4)-nh
            nk=NumGlu(5)-nj

            rIn = na+1
            rOut= NumGlu(1)+NumGlu(2)+ne
            Eps1(:) = cur_g_2f(Gluons(rIn:rOut),Quarks(1:2),(/1+nb+NumGlu(2)+ne,nb,NumGlu(2),ne/))
            PMom3 = SumMom(Gluons,rIn,rOut) + Quarks(1)%Mom + Quarks(2)%Mom
            PropFac3 = (0d0,-1d0)/sc_(PMom3,PMom3)! here was a bug: minus sign was missing
            if( abs(sc_(PMom3,PMom3)).lt.PropCut ) cycle
            Eps1 = Eps1*PropFac3
            
            rIn = NumGlu(1)+NumGlu(2)+ne+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+nh
            Sca2 = cur_s_2s(Gluons(rIn:rOut),Scalars(3:3),(/nf+nh,nf,nh/))
            PMom1 = SumMom(Gluons,rIn,rOut) + Scalars(3)%Mom
            if( nf.ge.1 .or. nh.ge.1 ) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1) - Scalars(3)%Mass2)
               if( abs(sc_(PMom1,PMom1) - Scalars(3)%Mass2).lt.PropCut ) cycle
               Sca2 = Sca2*PropFac1
            endif

            rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+nh+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+nj
            Sca1 = cur_s_2s(Gluons(rIn:rOut),Scalars(4:4),(/ni+nj,ni,nj/))
            PMom2 = SumMom(Gluons,rIn,rOut) + Scalars(4)%Mom
            if( ni.ge.1 .or. nj.ge.1 ) then
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Scalars(4)%Mass2)
               if( abs(sc_(PMom2,PMom2) - Scalars(4)%Mass2).lt.PropCut ) cycle
               Sca1 = Sca1*PropFac2
            endif

            Eps2(:) = vggss(Dv,Eps1) * Sca1*Sca2


            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpMom1(:) = PMom1(:)+PMom2(:)+PMom3(:)
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps2(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+nj+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps3(:) = cur_g(TmpGluons(1:counter-1),1+na+nk+1)

            if( na.ge.1 .or. nk.ge.1 ) then
               PropFac1 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
               if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
               Eps3 = Eps3*PropFac1
            endif

            Res = Res + Eps3

         enddo
         enddo
         enddo
         enddo



return
END FUNCTION





! this current only appears in AmpType 3
FUNCTION cur_g_sffs(Gluons,Scalars,Quarks,NumGlu) result(res)           ! Gluons(:) does not include the OFF-shell gluon,  however NumGlu is the number of all gluons
implicit none                                                           ! checked gauge inv. for one gluon   
integer,intent(in) :: NumGlu(0:5)
type(PtrToParticle) :: Gluons(1:),Scalars(1:2),Quarks(3:4)
integer :: na,nb,nc,nd,ne,nf,ng,nh,ni,nj,nk
integer :: rIn,rOut
integer :: counter,i
complex(8) :: res(Dv)
type(PtrToParticle) :: TmpGluons(1:2+NumGlu(1)+NumGlu(3)+NumGlu(5))
complex(8),target :: TmpMom1(1:Dv),TmpMom2(1:Dv)
integer,target :: TmpExtRef1,TmpExtRef2
complex(8),target :: Eps1(1:Dv)
complex(8),target :: Eps2(1:Dv)
complex(8) :: Eps3(1:Dv)
complex(8) :: Sca1,Sca2
complex(8) :: u1(1:Ds),ubar2(1:Ds)
complex(8) :: PropFac1,PropFac2,PropFac3,PropFac4
complex(8) :: PMom1(1:Dv)
complex(8) :: PMom2(1:Dv)
complex(8) :: PMom3(1:Dv)
complex(8) :: PMom4(1:Dv)

!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-1-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4)-NumGlu(5).ne.0 ) print *, "wrong number of gluons in cur_g_ssff"
!DEC$ ENDIF

      res = (0d0,0d0)
!        type (1)
         do na=0,NumGlu(1)
         do nc=0,NumGlu(2)
         do ne=0,NumGlu(5)
            nb=NumGlu(1)-na
            nd=NumGlu(2)-nc
            nf=NumGlu(5)-ne

            rIn = NumGlu(1)+nc+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne
            Sca1 = cur_s_sffs(Gluons(rIn:rOut),Scalars(2:2),Quarks(3:4),(/nd+NumGlu(3)+NumGlu(4)+ne,nd,NumGlu(3),NumGlu(4),ne/))
            PMom2  = SumMom(Gluons,rIn,rOut)  + Scalars(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom
            PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Scalars(1)%Mass2)
            if( abs(sc_(PMom2,PMom2) - Scalars(1)%Mass2).lt.PropCut ) cycle
            Sca1 = Sca1*PropFac2

            rIn = na+1
            rOut= NumGlu(1)+nc
            Sca2 = cur_s_2s(Gluons(rIn:rOut),Scalars(1:1),(/nb+nc,nb,nc/))
            PMom3 = SumMom(Gluons,rIn,rOut) + Scalars(1)%Mom
            if( nb.ge.1 .or. nc.ge.1 ) then
               PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Scalars(1)%Mass2)
               if( abs(sc_(PMom3,PMom3) - Scalars(1)%Mass2).lt.PropCut ) cycle
               Sca2 = Sca2*PropFac3
            endif

            if( Scalars(1)%PartType.gt.0 ) then
                Eps1(:) = vbss(Dv,PMom3(:),PMom2(:)) * Sca1*Sca2!  here was a bug: P2,P3
            else
                Eps1(:) = vsbs(Dv,PMom3(:),PMom2(:)) * Sca1*Sca2!  here was a bug: P2,P3
            endif
            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpMom1(:) = PMom2(:)+PMom3(:)
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps1(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps2(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+1)

            if( na.ge.1 .or. nf.ge.1 ) then
               PropFac1 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
               if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
               Eps2 = Eps2*PropFac1
            endif

            Res = Res + Eps2
         enddo
         enddo
         enddo


!        type (2)
         do na=0,NumGlu(1)
         do nc=0,NumGlu(4)
         do ne=0,NumGlu(5)
            nb=NumGlu(1)-na
            nd=NumGlu(4)-nc
            nf=NumGlu(5)-ne

            rIn = na+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+nc
            Sca1 = cur_s_ssff(Gluons(rIn:rOut),Scalars(1:1),Quarks(3:4),(/nb+NumGlu(2)+NumGlu(3)+nc,nb,NumGlu(2),NumGlu(3),nc/))
            PMom2 =  SumMom(Gluons,rIn,rOut) + Scalars(1)%Mom + Quarks(3)%Mom + Quarks(4)%Mom
            PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Scalars(2)%Mass2)
            if( abs(sc_(PMom2,PMom2) - Scalars(2)%Mass2).lt.PropCut ) cycle
            Sca1 = Sca1*PropFac2

            rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+nc+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne
            Sca2 = cur_s_2s(Gluons(rIn:rOut),Scalars(2:2),(/nd+ne,nd,ne/))
            PMom3  = SumMom(Gluons,rIn,rOut) + Scalars(2)%Mom
            if( nd.ge.1 .or. ne.ge.1 ) then
               PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Scalars(2)%Mass2)
               if( abs(sc_(PMom3,PMom3) - Scalars(2)%Mass2).lt.PropCut ) cycle
               Sca2 = Sca2*PropFac3
            endif

            if( Scalars(2)%PartType.gt.0 ) then
                Eps1(:) = vsbs(Dv,PMom2(:),PMom3(:)) * Sca1*Sca2
            else
                Eps1(:) = vbss(Dv,PMom2(:),PMom3(:)) * Sca1*Sca2
            endif


            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpMom1(:) = PMom2(:)+PMom3(:)
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps1(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps2(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+1)

            if( na.ge.1 .or. nf.ge.1 ) then
               PropFac1 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
               if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
               Eps2 = Eps2*PropFac1
            endif

            Res = Res + Eps2
         enddo
         enddo
         enddo




!        type (3)
         do na=0,NumGlu(1)
         do nc=0,NumGlu(2)
         do ne=0,NumGlu(4)
         do nj=0,NumGlu(5)
            nb=NumGlu(1)-na
            nd=NumGlu(2)-nc
            nf=NumGlu(4)-ne
            nk=NumGlu(5)-nj

            rIn = na+1
            rOut= NumGlu(1)+nc
            Sca2 = cur_s_2s(Gluons(rIn:rOut),Scalars(1:1),(/nb+nc,nb,nc/))
            PMom1 = SumMom(Gluons,rIn,rOut) + Scalars(1)%Mom
            if( nb.ge.1 .or. nc.ge.1 ) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1) - Scalars(1)%Mass2)
               if( abs(sc_(PMom1,PMom1) - Scalars(1)%Mass2).lt.PropCut ) cycle
               Sca2 = Sca2*PropFac1
            endif

            rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+ne+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+nj
            Sca1 = cur_s_2s(Gluons(rIn:rOut),Scalars(2:2),(/nf+nj,nf,nj/))
            PMom2 = SumMom(Gluons,rIn,rOut) + Scalars(2)%Mom
            if( nf.ge.1 .or. nj.ge.1 ) then
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Scalars(2)%Mass2)
               if( abs(sc_(PMom2,PMom2) - Scalars(2)%Mass2).lt.PropCut ) cycle
               Sca1 = Sca1*PropFac2
            endif

            rIn = NumGlu(1)+nc+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+ne
            Eps1(:) = cur_g_2f(Gluons(rIn:rOut),Quarks(3:4),(/1+nd+NumGlu(3)+ne,nd,NumGlu(3),ne/))
            PMom3 = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom + Quarks(4)%Mom
            PropFac3 = (0d0,-1d0)/sc_(PMom3,PMom3)! here was a bug: minus sign was missing
            if( abs(sc_(PMom3,PMom3)).lt.PropCut ) cycle
            Eps1 = Eps1*PropFac3


            Eps2(:) = vgsgs(Dv,Eps1) * Sca1*Sca2


            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpMom1(:) = PMom1(:)+PMom2(:)+PMom3(:)
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps2(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+nj+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps3(:) = cur_g(TmpGluons(1:counter-1),1+na+nk+1)

            if( na.ge.1 .or. nk.ge.1 ) then
               PropFac1 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
               if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
               Eps3 = Eps3*PropFac1
            endif

            Res = Res + Eps3
         enddo
         enddo
         enddo
         enddo

return
END FUNCTION








FUNCTION cur_s_ssffss(Gluons,Scalars,Quarks,NumGlu) result(res)
implicit none
integer :: NumGlu(0:6)
type(PtrToParticle) :: Gluons(1:),Quarks(3:4),Scalars(2:4)      ! off-shell scalar is not included
integer,target :: TmpExtRef
complex(8) :: res,tmp
complex(8),target :: Sca1(1:1)
complex(8) :: eps1(1:Dv)
complex(8) :: eps2(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(4)),TmpScalar(1:3)
complex(8) :: PropFac1,PropFac2
complex(8),target :: pmom1(1:Dv)
complex(8) :: pmom2(1:Dv)
integer :: n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,n5a,n5b,n6a,n6b
integer :: rIn,rOut,i,counter


!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4)-NumGlu(5)-NumGlu(6).ne.0 ) print *, "wrong number of gluons in cur_s_ssffss"
!DEC$ ENDIF


    res = (0d0,0d0)

!   (A)
      do n2a=0,NumGlu(2)
      do n6a=0,NumGlu(6)
         n2b = NumGlu(2)-n2a
         n6b = NumGlu(6)-n6a

         rIn =NumGlu(1)+n2a+1
         rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
         Eps2 = cur_g_ffss(Gluons(rIn:rOut),Scalars(3:4),Quarks(3:4),(/1+n2b+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a,n2b,NumGlu(3),NumGlu(4),NumGlu(5),n6a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom + Quarks(4)%Mom+ Scalars(3)%Mom + Scalars(4)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut) cycle
         Eps2 = Eps2*PropFac1
         do n1a=0,NumGlu(1)
            n1b = NumGlu(1)-n1a
            rIn =n1a+1
            rOut=NumGlu(1)+n2a
            sca1(1) = cur_s_2s(Gluons(rIn:rOut),Scalars(2:2),(/n1b+n2a,n1b,n2a/) )
            PMom2(:) = Scalars(2)%Mom + SumMom(Gluons,rIn,rOut)
            if(n1b.ge.1 .or. n2a.ge.1) then
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Scalars(2)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Scalars(2)%Mass2).lt.PropCut ) cycle
               sca1(1) = sca1(1)*PropFac2
            endif

           if( Scalars(2)%PartType.gt.0 ) then
              sca1(1) = csg(Eps2,PMom1,PMom2) * sca1(1)
           else
              sca1(1) = cbsg(Eps2,PMom1,PMom2) * sca1(1)
           endif

            rIn = n1a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            PMom1 = SumMom(Gluons,rIn,rOut) + Scalars(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom + Scalars(3)%Mom + Scalars(4)%Mom
            if(n1a.ge.1 .or. n6b.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Scalars(2)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Scalars(2)%Mass2).lt.PropCut ) cycle
               sca1(1) = sca1(1)*PropFac2
            endif

            TmpScalar(1)%Mom  => PMom1(:)
            TmpScalar(1)%Pol  => sca1
            TmpScalar(1)%Mass => Scalars(2)%Mass
            TmpScalar(1)%Mass2=> Scalars(2)%Mass2
            TmpExtRef = -1
            TmpScalar(1)%ExtRef => TmpExtRef
            TmpScalar(1)%PartType => Scalars(2)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp = cur_s_2s(TmpGluons(1:counter-1),TmpScalar(1:1),(/counter-1,n1a,n6b/) )

            Res = Res + tmp
         enddo
      enddo
      enddo





!   (B)
      do n1a=0,NumGlu(1)
      do n4a=0,NumGlu(4)
         n1b = NumGlu(1)-n1a
         n4b = NumGlu(4)-n4a

            rIn = n1a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
            sca1 = cur_s_ssff(Gluons(rIn:rOut),Scalars(2:2),Quarks(3:4),(/n1b+NumGlu(2)+NumGlu(3)+n4a,n1b,NumGlu(2),NumGlu(3),n4a/))
            PMom1  = SumMom(Gluons,rIn,rOut)  + Scalars(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom
            PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1) - Scalars(2)%Mass2)
            if( abs(sc_(PMom1,PMom1) - Scalars(2)%Mass2).lt.PropCut ) cycle
            sca1 = sca1 * PropFac1

            TmpScalar(1)%Mom  => PMom1(:)
            TmpScalar(1)%Pol  => sca1
            TmpScalar(1)%Mass => Scalars(2)%Mass
            TmpScalar(1)%Mass2=> Scalars(2)%Mass2
            TmpExtRef = -1
            TmpScalar(1)%ExtRef => TmpExtRef
            TmpScalar(1)%PartType => Scalars(2)%PartType        
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
            rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            call CopyParticlePtr(Scalars(3),TmpScalar(2))
            call CopyParticlePtr(Scalars(4),TmpScalar(3))
            tmp = cur_s_4s(TmpGluons(1:counter-1),TmpScalar(1:3),(/counter-1,n1a,n4b,NumGlu(5),NumGlu(6)/) )

            Res = Res + tmp
!             Res3(:) = Res3(:) + tmp(:)
      enddo
      enddo







!   (C)
      do n1a=0,NumGlu(1)
      do n2a=0,NumGlu(2)
      do n6a=0,NumGlu(6)
         n1b = NumGlu(1)-n1a
         n2b = NumGlu(2)-n2a
         n6b = NumGlu(6)-n6a

            rIn =n1a+1
            rOut=NumGlu(1)+n2a
            sca1 = cur_s_2s(Gluons(rIn:rOut),Scalars(2:2),(/n1b+n2a,n1b,n2a/) )
            PMom2(:) = Scalars(2)%Mom + SumMom(Gluons,rIn,rOut)
            if(n1b.ge.1 .or. n2a.ge.1) then
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Scalars(2)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Scalars(2)%Mass2).lt.PropCut ) cycle
               sca1 = sca1*PropFac2
            endif

            do n4a=0,NumGlu(4)
                n4b = NumGlu(4)-n4a

                rIn = NumGlu(1)+n2a+1
                rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
                Eps1 = cur_g_2f(Gluons(rIn:rOut),Quarks(3:4),(/1+n2b+NumGlu(3)+n4a,n2b,NumGlu(3),n4a/))
                PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom + Quarks(4)%Mom
                PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
                if( abs(sc_(PMom1,PMom1)).lt.PropCut) cycle
                Eps1 = Eps1*PropFac1

                rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                Eps2 = cur_g_2s(Gluons(rIn:rOut),Scalars(3:4),(/1+n4b+NumGlu(5)+n6a,n4b,NumGlu(5),n6a/))
                PMom2(:) = SumMom(Gluons,rIn,rOut) + Scalars(3)%Mom + Scalars(4)%Mom
                PropFac2 = (0d0,-1d0)/sc_(PMom2,PMom2)
                if( abs(sc_(PMom2,PMom2)).lt.PropCut) cycle
                Eps2 = Eps2*PropFac2

                sca1 = csgg(Eps1,Eps2) * sca1

                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                  PMom1 = SumMom(Gluons,rIn,rOut) + Scalars(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom + Scalars(3)%Mom + Scalars(4)%Mom
                  if(n1a.ge.1 .or. n6b.ge.1) then
                    PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Scalars(2)%Mass2)
                    if( abs(sc_(PMom1,PMom1)-Scalars(2)%Mass2).lt.PropCut ) cycle
                    sca1 = sca1*PropFac2
                  endif

                  TmpScalar(1)%Mom  => PMom1(:)
                  TmpScalar(1)%Pol  => sca1
                  TmpScalar(1)%Mass => Scalars(2)%Mass
                  TmpScalar(1)%Mass2=> Scalars(2)%Mass2
                  TmpExtRef = -1
                  TmpScalar(1)%ExtRef => TmpExtRef
                  TmpScalar(1)%PartType => Scalars(2)%PartType
                  counter=1
                  rIn =1
                  rOut=n1a
                  do i=rIn,rOut
                    call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                    counter=counter+1
                  enddo
                  rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
                  rOut=NumGlu(0)
                  do i=rIn,rOut
                    call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                    counter=counter+1
                  enddo
                  tmp = cur_s_2s(TmpGluons(1:counter-1),TmpScalar(1:1),(/counter-1,n1a,n6b/) )

                  Res = Res + tmp
      !             Res1(:) = Res1(:) + tmp(:)
      !             print *, "1",tmp(:)
            enddo
      enddo
      enddo
      enddo



return
END FUNCTION






! this current only appears in AmpType 3
FUNCTION cur_s_sffsss(Gluons,Scalars,Quarks,NumGlu) result(res)
implicit none
integer :: NumGlu(0:6)
type(PtrToParticle) :: Gluons(1:),Quarks(3:4),Scalars(2:4)      ! off-shell scalar is not included
integer,target :: TmpExtRef
complex(8) :: res,tmp
complex(8),target :: Sca1(1:1)
complex(8) :: eps1(1:Dv)
complex(8) :: eps2(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(4)),TmpScalar(1:3)
complex(8) :: PropFac1,PropFac2
complex(8),target :: pmom1(1:Dv)
complex(8) :: pmom2(1:Dv)
integer :: n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,n5a,n5b,n6a,n6b
integer :: rIn,rOut,i,counter


!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4)-NumGlu(5)-NumGlu(6).ne.0 ) print *, "wrong number of gluons in cur_s_sffsss"
!DEC$ ENDIF


if( NumGlu(0).ne.0 ) call Error("cur_s_sffsss works only with no gluons so far")

    res = (0d0,0d0)

!   type (A)
      do n3a=0,NumGlu(3)
      do n6a=0,NumGlu(6)
         n3b = NumGlu(3)-n3a
         n6b = NumGlu(6)-n6a

            rIn = NumGlu(1)+NumGlu(2)+n3a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6b
            sca1 = cur_s_sffs(Gluons(rIn:rOut),Scalars(2:2),Quarks(3:4),(/n3b+NumGlu(4)+NumGlu(5)+n6a, n3b,NumGlu(4),NumGlu(5),n6a/))
            PMom1  = SumMom(Gluons,rIn,rOut)  + Scalars(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom
            PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1) - Scalars(2)%Mass2)
            if( abs(sc_(PMom1,PMom1) - Scalars(2)%Mass2).lt.PropCut ) cycle
            sca1 = sca1 * PropFac1

            TmpScalar(1)%Mom  => PMom1(:)
            TmpScalar(1)%Pol  => sca1
            TmpScalar(1)%Mass => Scalars(2)%Mass
            TmpScalar(1)%Mass2=> Scalars(2)%Mass2
            TmpExtRef = -1
            TmpScalar(1)%ExtRef => TmpExtRef
            TmpScalar(1)%PartType => Scalars(2)%PartType      
            counter=1
            rIn =1
            rOut=NumGlu(1)+NumGlu(2)+n3a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+NumGlu(6)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp = cur_s_4s(TmpGluons(1:counter-1),(/TmpScalar(1),Scalars(3),Scalars(4)/),(/NumGlu(1)+NumGlu(2)+n3a+n6b,NumGlu(1),NumGlu(2),n3a,n6b/) )
            Res = Res + tmp
      enddo
      enddo




!   type (B)
      do n3a=0,NumGlu(3)
      do n6a=0,NumGlu(6)
         n3b = NumGlu(3)-n3a
         n6b = NumGlu(6)-n6a

            rIn = NumGlu(1)+NumGlu(2)+n3a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6b
            sca1 = cur_s_4s(Gluons(rIn:rOut),Scalars(2:4),(/n3b+NumGlu(4)+NumGlu(5)+n6a, n3b,NumGlu(4),NumGlu(5),n6a/))
            PMom1  = SumMom(Gluons,rIn,rOut)  + Scalars(2)%Mom + Scalars(3)%Mom + Scalars(4)%Mom
            PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1) - Scalars(2)%Mass2)
            if( abs(sc_(PMom1,PMom1) - Scalars(2)%Mass2).lt.PropCut ) cycle
            sca1 = sca1 * PropFac1

            TmpScalar(1)%Mom  => PMom1(:)
            TmpScalar(1)%Pol  => sca1
            TmpScalar(1)%Mass => Scalars(2)%Mass
            TmpScalar(1)%Mass2=> Scalars(2)%Mass2
            TmpExtRef = -1
            TmpScalar(1)%ExtRef => TmpExtRef
            TmpScalar(1)%PartType => Scalars(2)%PartType        
            counter=1
            rIn =1
            rOut=NumGlu(1)+NumGlu(2)+n3a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+NumGlu(6)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp = cur_s_sffs(TmpGluons(1:counter-1),TmpScalar(1),Quarks(3:4),(/NumGlu(1)+NumGlu(2)+n3a+n6b,NumGlu(1),NumGlu(2),n3a,n6b/) )
            Res = Res + tmp
      enddo
      enddo







!   type (C)
      do n1a=0,NumGlu(1)
      do n5a=0,NumGlu(5)
      do n6a=0,NumGlu(6)
         n1b = NumGlu(1)-n1a
         n5b = NumGlu(5)-n5a
         n6b = NumGlu(6)-n6a

            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a+1
            rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            sca1 = cur_s_2s(Gluons(rIn:rOut),Scalars(2:2),(/n5b+n6a,n5b,n6a/) )! here was a bug n5b,n6a

            PMom2(:) = Scalars(2)%Mom + SumMom(Gluons,rIn,rOut)
            if(n5b.ge.1 .or. n6a.ge.1) then
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Scalars(2)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Scalars(2)%Mass2).lt.PropCut ) cycle
               sca1 = sca1*PropFac2
            endif

            do n3a=0,NumGlu(3)
                n3b = NumGlu(3)-n3a

                rIn = n1a+1
                rOut= NumGlu(1)+NumGlu(2)+n3a
                Eps1 = cur_g_2f(Gluons(rIn:rOut),Quarks(3:4),(/1+n1b+NumGlu(2)+n3a,n1b,NumGlu(2),n3a/))
                PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom + Quarks(4)%Mom
                PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
                if( abs(sc_(PMom1,PMom1)).lt.PropCut) cycle
                Eps1 = Eps1*PropFac1

                rIn = NumGlu(1)+NumGlu(2)+n3a+1
                rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a
                Eps2 = cur_g_2s(Gluons(rIn:rOut),Scalars(3:4),(/1+n3b+NumGlu(4)+n5a,n3b,NumGlu(4),n5a/))
                PMom2(:) = SumMom(Gluons,rIn,rOut) + Scalars(3)%Mom + Scalars(4)%Mom
                PropFac2 = (0d0,-1d0)/sc_(PMom2,PMom2)
                if( abs(sc_(PMom2,PMom2)).lt.PropCut) cycle
                Eps2 = Eps2*PropFac2

                sca1 = cgsg(Eps1,Eps2) * sca1

                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                  PMom1 = SumMom(Gluons,rIn,rOut) + Scalars(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom + Scalars(3)%Mom + Scalars(4)%Mom
                  if(n1a.ge.1 .or. n6b.ge.1) then
                    PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Scalars(2)%Mass2)
                    if( abs(sc_(PMom1,PMom1)-Scalars(2)%Mass2).lt.PropCut ) cycle
                    sca1 = sca1*PropFac2
                  endif

                  TmpScalar(1)%Mom  => PMom1(:)
                  TmpScalar(1)%Pol  => sca1
                  TmpScalar(1)%Mass => Scalars(2)%Mass
                  TmpScalar(1)%Mass2=> Scalars(2)%Mass2
                  TmpExtRef = -1
                  TmpScalar(1)%ExtRef => TmpExtRef
                  TmpScalar(1)%PartType => Scalars(2)%PartType
                  counter=1
                  rIn =1
                  rOut=n1a
                  do i=rIn,rOut
                    call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                    counter=counter+1
                  enddo
                  rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
                  rOut=NumGlu(0)
                  do i=rIn,rOut
                    call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                    counter=counter+1
                  enddo
                  tmp = cur_s_2s(TmpGluons(1:counter-1),TmpScalar(1:1),(/counter-1,n1a,n6b/) )

                  Res = Res + tmp
            enddo
      enddo
      enddo
      enddo


return
END FUNCTION










! this current only appears in scalar closed loop
FUNCTION cur_s_sffsss_CLOSEDLOOPCONTRIB(Gluons,Scalars,Quarks,NumGlu) result(res)
implicit none
integer :: NumGlu(0:6)
type(PtrToParticle) :: Gluons(1:),Quarks(3:4),Scalars(2:4)      ! off-shell scalar is not included
integer,target :: TmpExtRef
complex(8) :: res,tmp
complex(8),target :: Sca1(1:1)
complex(8) :: eps1(1:Dv)
complex(8) :: eps2(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(4)),TmpScalar(1:3)
complex(8) :: PropFac1,PropFac2
complex(8),target :: pmom1(1:Dv)
complex(8) :: pmom2(1:Dv)
integer :: n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,n5a,n5b,n6a,n6b
integer :: rIn,rOut,i,counter


!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4)-NumGlu(5)-NumGlu(6).ne.0 ) print *, "wrong number of gluons in cur_s_sffsss_CLOSEDLOOPCONTRIB"
!DEC$ ENDIF


    res = (0d0,0d0)

!   type (A)
      do n1a=0,NumGlu(1)
      do n5a=0,NumGlu(5)
         n1b = NumGlu(1)-n1a
         n5b = NumGlu(5)-n5a

         rIn =n1a+1
         rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a
         Eps2 = cur_g_ffss(Gluons(rIn:rOut),Scalars(2:3),Quarks(3:4),(/1+n1b+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a, n1b,NumGlu(2),NumGlu(3),NumGlu(4),n5a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom + Quarks(4)%Mom+ Scalars(2)%Mom + Scalars(3)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut) cycle
         Eps2 = Eps2*PropFac1
         do n6a=0,NumGlu(6)
            n6b = NumGlu(6)-n6a
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a+1
            rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            sca1(1) = cur_s_2s(Gluons(rIn:rOut),Scalars(4:4),(/n5b+n6a,n5b,n6a/) )
            PMom2(:) = Scalars(4)%Mom + SumMom(Gluons,rIn,rOut)
            if(n5b.ge.1 .or. n6a.ge.1) then
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Scalars(4)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Scalars(4)%Mass2).lt.PropCut ) cycle
               sca1(1) = sca1(1)*PropFac2
            endif

            if( Scalars(4)%PartType.gt.0 ) then
                sca1(1) = cgs(Eps2,PMom1,PMom2) * sca1(1)
            else
                sca1(1) = cgbs(Eps2,PMom1,PMom2) * sca1(1)
            endif

            rIn = n1a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            PMom1 = SumMom(Gluons,rIn,rOut) + Scalars(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom + Scalars(3)%Mom + Scalars(4)%Mom
            if(n1a.ge.1 .or. n6b.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Scalars(4)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Scalars(4)%Mass2).lt.PropCut ) cycle
               sca1(1) = sca1(1)*PropFac2
            endif

            TmpScalar(1)%Mom  => PMom1(:)
            TmpScalar(1)%Pol  => sca1
            TmpScalar(1)%Mass => Scalars(4)%Mass
            TmpScalar(1)%Mass2=> Scalars(4)%Mass2
            TmpExtRef = -1
            TmpScalar(1)%ExtRef => TmpExtRef
            TmpScalar(1)%PartType => Scalars(4)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp = cur_s_2s(TmpGluons(1:counter-1),TmpScalar(1:1),(/counter-1,n1a,n6b/) )

            Res = Res + tmp
         enddo
      enddo
      enddo





!   type (B)
      do n3a=0,NumGlu(3)
      do n6a=0,NumGlu(6)
         n3b = NumGlu(3)-n3a
         n6b = NumGlu(6)-n6a

            rIn = NumGlu(1)+NumGlu(2)+n3a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6b
            sca1 = cur_s_4s(Gluons(rIn:rOut),Scalars(2:4),(/n3b+NumGlu(4)+NumGlu(5)+n6a, n3b,NumGlu(4),NumGlu(5),n6a/),tag_f=3)
            PMom1  = SumMom(Gluons,rIn,rOut)  + Scalars(2)%Mom + Scalars(3)%Mom + Scalars(4)%Mom
            PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1) - Scalars(4)%Mass2)
            if( abs(sc_(PMom1,PMom1) - Scalars(4)%Mass2).lt.PropCut ) cycle
            sca1 = sca1 * PropFac1

            TmpScalar(1)%Mom  => PMom1(:)
            TmpScalar(1)%Pol  => sca1
            TmpScalar(1)%Mass => Scalars(4)%Mass
            TmpScalar(1)%Mass2=> Scalars(4)%Mass2
            TmpExtRef = -1
            TmpScalar(1)%ExtRef => TmpExtRef
            TmpScalar(1)%PartType => Scalars(4)%PartType        
            counter=1
            rIn =1
            rOut=NumGlu(1)+NumGlu(2)+n3a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+NumGlu(6)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp = cur_s_sffs(TmpGluons(1:counter-1),TmpScalar(1),Quarks(3:4),(/NumGlu(1)+NumGlu(2)+n3a+n6b,NumGlu(1),NumGlu(2),n3a,n6b/) )
            Res = Res + tmp
      enddo
      enddo







!   type (C)
      do n1a=0,NumGlu(1)
      do n5a=0,NumGlu(5)
      do n6a=0,NumGlu(6)
         n1b = NumGlu(1)-n1a
         n5b = NumGlu(5)-n5a
         n6b = NumGlu(6)-n6a

            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a+1
            rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            sca1 = cur_s_2s(Gluons(rIn:rOut),Scalars(4:4),(/n5b+n6a,n5b,n6a/) )! here was a bug n5b,n6a, and another bug: Scalars(2-->4)

            PMom2(:) = Scalars(4)%Mom + SumMom(Gluons,rIn,rOut)
            if(n5b.ge.1 .or. n6a.ge.1) then
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Scalars(4)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Scalars(4)%Mass2).lt.PropCut ) cycle
               sca1 = sca1*PropFac2
            endif

            do n3a=0,NumGlu(3)
                n3b = NumGlu(3)-n3a

                rIn = n1a+1
                rOut= NumGlu(1)+NumGlu(2)+n3a
                Eps1 = cur_g_2f(Gluons(rIn:rOut),Quarks(3:4),(/1+n1b+NumGlu(2)+n3a,n1b,NumGlu(2),n3a/))
                PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom + Quarks(4)%Mom
                PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
                if( abs(sc_(PMom1,PMom1)).lt.PropCut) cycle
                Eps1 = Eps1*PropFac1

                rIn = NumGlu(1)+NumGlu(2)+n3a+1
                rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a
                Eps2 = cur_g_2s(Gluons(rIn:rOut),Scalars(2:3),(/1+n3b+NumGlu(4)+n5a,n3b,NumGlu(4),n5a/))
                PMom2(:) = SumMom(Gluons,rIn,rOut) + Scalars(2)%Mom + Scalars(3)%Mom
                PropFac2 = (0d0,-1d0)/sc_(PMom2,PMom2)
                if( abs(sc_(PMom2,PMom2)).lt.PropCut) cycle
                Eps2 = Eps2*PropFac2

                sca1 = cggs(Eps1,Eps2) * sca1

                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                  PMom1 = SumMom(Gluons,rIn,rOut) + Scalars(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom + Scalars(3)%Mom + Scalars(4)%Mom
                  if(n1a.ge.1 .or. n6b.ge.1) then
                    PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Scalars(2)%Mass2)
                    if( abs(sc_(PMom1,PMom1)-Scalars(2)%Mass2).lt.PropCut ) cycle
                    sca1 = sca1*PropFac2
                  endif

                  TmpScalar(1)%Mom  => PMom1(:)
                  TmpScalar(1)%Pol  => sca1
                  TmpScalar(1)%Mass => Scalars(4)%Mass
                  TmpScalar(1)%Mass2=> Scalars(4)%Mass2
                  TmpExtRef = -1
                  TmpScalar(1)%ExtRef => TmpExtRef
                  TmpScalar(1)%PartType => Scalars(4)%PartType
                  counter=1
                  rIn =1
                  rOut=n1a
                  do i=rIn,rOut
                    call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                    counter=counter+1
                  enddo
                  rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
                  rOut=NumGlu(0)
                  do i=rIn,rOut
                    call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                    counter=counter+1
                  enddo
                  tmp = cur_s_2s(TmpGluons(1:counter-1),TmpScalar(1:1),(/counter-1,n1a,n6b/) )

                  Res = Res + tmp
            enddo
      enddo
      enddo
      enddo


return
END FUNCTION



! this current only appears in closed scalar loop
FUNCTION cur_s_sssffs_CLOSEDLOOPCONTRIB(Gluons,Scalars,Quarks,NumGlu) result(res)
implicit none
integer :: NumGlu(0:6)
type(PtrToParticle) :: Gluons(1:),Quarks(3:4),Scalars(2:4)      ! off-shell scalar is not included
integer,target :: TmpExtRef
complex(8) :: res,tmp
complex(8),target :: Sca1(1:1)
complex(8) :: eps1(1:Dv)
complex(8) :: eps2(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(4)),TmpScalar(1:3)
complex(8) :: PropFac1,PropFac2
complex(8),target :: pmom1(1:Dv)
complex(8) :: pmom2(1:Dv)
integer :: n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,n5a,n5b,n6a,n6b
integer :: rIn,rOut,i,counter


!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4)-NumGlu(5)-NumGlu(6).ne.0 ) print *, "wrong number of gluons in cur_s_sssffs_CLOSEDLOOPCONTRIB"
!DEC$ ENDIF


    res = (0d0,0d0)

!   type (A)
      do n1a=0,NumGlu(1)
      do n5a=0,NumGlu(5)
         n1b = NumGlu(1)-n1a
         n5b = NumGlu(5)-n5a

         rIn =n1a+1
         rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a
         Eps2 = cur_g_ssff(Gluons(rIn:rOut),Scalars(2:3),Quarks(3:4),(/1+n1b+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a, n1b,NumGlu(2),NumGlu(3),NumGlu(4),n5a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom + Quarks(4)%Mom+ Scalars(2)%Mom + Scalars(3)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut) cycle
         Eps2 = Eps2*PropFac1
         do n6a=0,NumGlu(6)
            n6b = NumGlu(6)-n6a
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a+1
            rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            sca1(1) = cur_s_2s(Gluons(rIn:rOut),Scalars(4:4),(/n5b+n6a,n5b,n6a/) )
            PMom2(:) = Scalars(4)%Mom + SumMom(Gluons,rIn,rOut)
            if(n5b.ge.1 .or. n6a.ge.1) then
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Scalars(4)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Scalars(4)%Mass2).lt.PropCut ) cycle
               sca1(1) = sca1(1)*PropFac2
            endif

            if( Scalars(4)%PartType.gt.0 ) then
                sca1(1) = cgs(Eps2,PMom1,PMom2) * sca1(1)
            else
                sca1(1) = cgbs(Eps2,PMom1,PMom2) * sca1(1)
            endif
            rIn = n1a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            PMom1 = SumMom(Gluons,rIn,rOut) + Scalars(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom + Scalars(3)%Mom + Scalars(4)%Mom
            if(n1a.ge.1 .or. n6b.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Scalars(4)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Scalars(4)%Mass2).lt.PropCut ) cycle
               sca1(1) = sca1(1)*PropFac2
            endif

            TmpScalar(1)%Mom  => PMom1(:)
            TmpScalar(1)%Pol  => sca1
            TmpScalar(1)%Mass => Scalars(4)%Mass
            TmpScalar(1)%Mass2=> Scalars(4)%Mass2
            TmpExtRef = -1
            TmpScalar(1)%ExtRef => TmpExtRef
            TmpScalar(1)%PartType => Scalars(4)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp = cur_s_2s(TmpGluons(1:counter-1),TmpScalar(1:1),(/counter-1,n1a,n6b/) )

            Res = Res + tmp
         enddo
      enddo
      enddo



!   type (B)
      do n3a=0,NumGlu(3)
      do n6a=0,NumGlu(6)
         n3b = NumGlu(3)-n3a
         n6b = NumGlu(6)-n6a

            rIn = NumGlu(1)+NumGlu(2)+n3a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6b
            sca1 = cur_s_sffs(Gluons(rIn:rOut),Scalars(4:4),Quarks(3:4),(/n3b+NumGlu(4)+NumGlu(5)+n6a, n3b,NumGlu(4),NumGlu(5),n6a/))
            PMom1  = SumMom(Gluons,rIn,rOut)  + Scalars(4)%Mom + Quarks(3)%Mom + Quarks(4)%Mom
            PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1) - Scalars(4)%Mass2)
            if( abs(sc_(PMom1,PMom1) - Scalars(4)%Mass2).lt.PropCut ) cycle
            sca1 = sca1 * PropFac1

            TmpScalar(1)%Mom  => PMom1(:)
            TmpScalar(1)%Pol  => sca1
            TmpScalar(1)%Mass => Scalars(4)%Mass
            TmpScalar(1)%Mass2=> Scalars(4)%Mass2
            TmpExtRef = -1
            TmpScalar(1)%ExtRef => TmpExtRef
            TmpScalar(1)%PartType => Scalars(4)%PartType        
            counter=1
            rIn =1
            rOut=NumGlu(1)+NumGlu(2)+n3a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+NumGlu(6)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp = cur_s_4s(TmpGluons(1:counter-1),(/Scalars(2),Scalars(3),TmpScalar(1)/),(/NumGlu(1)+NumGlu(2)+n3a+n6b,NumGlu(1),NumGlu(2),n3a,n6b/),tag_f=3)
            Res = Res + tmp
      enddo
      enddo






!   type (C)
      do n1a=0,NumGlu(1)
      do n5a=0,NumGlu(5)
      do n6a=0,NumGlu(6)
         n1b = NumGlu(1)-n1a
         n5b = NumGlu(5)-n5a
         n6b = NumGlu(6)-n6a

! ordering of cur_g_2f and cur_g_2s need to be exchanged if NumGlu >0

            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a+1
            rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            sca1 = cur_s_2s(Gluons(rIn:rOut),Scalars(4:4),(/n5b+n6a,n5b,n6a/) )! here was a bug n5b,n6a, and another bug: Scalars(2-->4)

            PMom2(:) = Scalars(4)%Mom + SumMom(Gluons,rIn,rOut)
            if(n5b.ge.1 .or. n6a.ge.1) then
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Scalars(4)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Scalars(4)%Mass2).lt.PropCut ) cycle
               sca1 = sca1*PropFac2
            endif

            do n3a=0,NumGlu(3)
                n3b = NumGlu(3)-n3a

                rIn = n1a+1
                rOut= NumGlu(1)+NumGlu(2)+n3a
                Eps1 = cur_g_2f(Gluons(rIn:rOut),Quarks(3:4),(/1+n1b+NumGlu(2)+n3a,n1b,NumGlu(2),n3a/))
                PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom + Quarks(4)%Mom
                PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
                if( abs(sc_(PMom1,PMom1)).lt.PropCut) cycle
                Eps1 = Eps1*PropFac1

                rIn = NumGlu(1)+NumGlu(2)+n3a+1
                rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a
                Eps2 = cur_g_2s(Gluons(rIn:rOut),Scalars(2:3),(/1+n3b+NumGlu(4)+n5a,n3b,NumGlu(4),n5a/))
                PMom2(:) = SumMom(Gluons,rIn,rOut) + Scalars(2)%Mom + Scalars(3)%Mom
                PropFac2 = (0d0,-1d0)/sc_(PMom2,PMom2)
                if( abs(sc_(PMom2,PMom2)).lt.PropCut) cycle
                Eps2 = Eps2*PropFac2

                sca1 = cggs(Eps1,Eps2) * sca1

                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                  PMom1 = SumMom(Gluons,rIn,rOut) + Scalars(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom + Scalars(3)%Mom + Scalars(4)%Mom
                  if(n1a.ge.1 .or. n6b.ge.1) then
                    PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Scalars(2)%Mass2)
                    if( abs(sc_(PMom1,PMom1)-Scalars(2)%Mass2).lt.PropCut ) cycle
                    sca1 = sca1*PropFac2
                  endif

                  TmpScalar(1)%Mom  => PMom1(:)
                  TmpScalar(1)%Pol  => sca1
                  TmpScalar(1)%Mass => Scalars(4)%Mass
                  TmpScalar(1)%Mass2=> Scalars(4)%Mass2
                  TmpExtRef = -1
                  TmpScalar(1)%ExtRef => TmpExtRef
                  TmpScalar(1)%PartType => Scalars(4)%PartType
                  counter=1
                  rIn =1
                  rOut=n1a
                  do i=rIn,rOut
                    call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                    counter=counter+1
                  enddo
                  rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
                  rOut=NumGlu(0)
                  do i=rIn,rOut
                    call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                    counter=counter+1
                  enddo
                  tmp = cur_s_2s(TmpGluons(1:counter-1),TmpScalar(1:1),(/counter-1,n1a,n6b/) )

                  Res = Res + tmp
            enddo
      enddo
      enddo
      enddo


return
END FUNCTION



! this current only appears in AmpType 3
FUNCTION cur_s_sssffs(Gluons,Scalars,Quarks,NumGlu) result(res)
implicit none
integer :: NumGlu(0:6)
type(PtrToParticle) :: Gluons(1:),Quarks(4:5),Scalars(2:4)      ! off-shell scalar is not included
integer,target :: TmpExtRef
complex(8) :: res,tmp
complex(8),target :: Sca0(1:1)
complex(8) :: Sca1
complex(8) :: EpsX(1:Dv)
complex(8) :: Eps2(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(6)),TmpScalar(1:1)
complex(8) :: PropFac1,PropFac2,PropFac4
complex(8),target :: pmom1(1:Dv)
complex(8) :: pmom2(1:Dv),pmom4(1:Dv)
integer :: n1a,n1b,n2a,n2b,n3a,n3b,n6a,n6b
integer :: rIn,rOut,i,counter


!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4)-NumGlu(5)-NumGlu(6).ne.0 ) print *, "wrong number of gluons in cur_s_sssffs"
!DEC$ ENDIF


      Res=(0d0,0d0)
      do n1a=0,NumGlu(1)
      do n2a=0,NumGlu(2)
      do n6a=0,NumGlu(6)
         n1b = NumGlu(1)-n1a
         n2b = NumGlu(2)-n2a
         n6b = NumGlu(6)-n6a

         rIn =NumGlu(1)+n2a+1
         rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
         Eps2 = cur_g_sffs(Gluons(rIn:rOut),Scalars(3:4),Quarks(4:5),(/1+n2b+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a,n2b,NumGlu(3),NumGlu(4),NumGlu(5),n6a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Scalars(3)%Mom + Scalars(4)%Mom + Quarks(4)%Mom + Quarks(5)%Mom

         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
         Eps2 = Eps2*PropFac1

            rIn =n1a+1
            rOut=NumGlu(1)+n2a
            Sca1 = cur_s_2s(Gluons(rIn:rOut),Scalars(2:2),(/n2a+n1b,n1b,n2a/))
            PMom2(:) = Scalars(2)%Mom + SumMom(Gluons,rIn,rOut)
            if(n1b.ge.1 .or. n2a.ge.1) then
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Scalars(2)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Scalars(2)%Mass2).lt.PropCut ) cycle
               Sca1 = Sca1*PropFac2
            endif

            if( Scalars(2)%PartType.gt.0 ) then
               Sca0 = csg(Eps2,PMom1,PMom2) * Sca1
            else
               Sca0 = cbsg(Eps2,PMom1,PMom2) * Sca1
            endif

            PMom1 = Scalars(2)%Mom + Scalars(3)%Mom + Scalars(4)%Mom + Quarks(4)%Mom + Quarks(5)%Mom + SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a)
            if(n1a.ge.1 .or. n6b.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Scalars(2)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Scalars(2)%Mass2).lt.PropCut ) cycle
               Sca0 = Sca0*PropFac1
            endif

            TmpScalar(1)%Mom  => PMom1(:)
            TmpScalar(1)%Pol  => Sca0
            TmpScalar(1)%Mass => Scalars(2)%Mass
            TmpScalar(1)%Mass2=> Scalars(2)%Mass2
            TmpExtRef = -1
            TmpScalar(1)%ExtRef => TmpExtRef
            TmpScalar(1)%PartType => Scalars(2)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp = cur_s_2s(TmpGluons(1:counter-1),TmpScalar(1:1),(/counter-1,n1a,n6b/) )
            Res = Res + tmp
      enddo
      enddo
      enddo



return
END FUNCTION








FUNCTION cur_s_4s(Gluons,Scalars,NumGlu,tag_f) result(res)!  checked gauge invariance for 2 gluons
implicit none
integer :: NumGlu(0:4)
type(PtrToParticle) :: Gluons(1:),Scalars(2:4)      ! off-shell scalar is not included
integer,target :: TmpExtRef
integer,optional :: tag_f
complex(8) :: res,tmp
complex(8),target :: Sca0(1:1)
complex(8) :: Sca1,Sca2,Sca3
complex(8) :: EpsX(1:Dv)
complex(8) :: Eps2(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(4)),TmpScalar(1:1)
complex(8) :: PropFac1,PropFac2,PropFac3,PropFac4
complex(8),target :: PMom1(1:Dv)
complex(8) :: PMom2(1:Dv),PMom3(1:Dv),PMom4(1:Dv)
integer :: n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,n1c,n2c,n4c,n3c
integer :: rIn,rOut,i,counter


!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4).ne.0 ) print *, "wrong number of gluons in cur_s_4s"
!DEC$ ENDIF


if( present(tag_f) ) then
if( tag_f.eq.3 .and. NumGlu(2)+NumGlu(4).eq.0 ) then! this calculates the inverted diagram with s1=s4 and s2=s3


! !     OLD CODE
!       Res=(0d0,0d0)
!          rIn =1
!          rOut=0
!          Eps2 = cur_g_2s(Gluons(rIn:rOut),Scalars(2:3),(/1,0,0,0/))
!          PMom1(:) = SumMom(Gluons,rIn,rOut) + Scalars(2)%Mom + Scalars(3)%Mom
!          PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
!          if( abs(sc_(PMom1,PMom1)).lt.PropCut ) return
!          Eps2 = Eps2*PropFac1
! 
!             rIn =1
!             rOut=0
!             Sca1 = cur_s_2s(Gluons(rIn:rOut),Scalars(4:4),(/0,0,0/))
!             PMom2(:) = Scalars(4)%Mom
!                   if( Scalars(4)%PartType.gt.0 ) then
!                       Sca0 = cgs(Eps2,PMom1,PMom2) * Sca1
!                   else
!                       Sca0 = cgbs(Eps2,PMom1,PMom2) * Sca1
!                   endif
! 
!             PMom1 = Scalars(2)%Mom+Scalars(3)%Mom+Scalars(4)%Mom
! 
!             TmpScalar(1)%Mom  => PMom1(:)
!             TmpScalar(1)%Pol  => Sca0
!             TmpScalar(1)%Mass => Scalars(4)%Mass
!             TmpScalar(1)%Mass2=> Scalars(4)%Mass2
!             TmpExtRef = -1
!             TmpScalar(1)%ExtRef => TmpExtRef
!             TmpScalar(1)%PartType => Scalars(4)%PartType
!             counter=1
!             rIn =1
!             rOut=0
!             do i=rIn,rOut
!               call CopyParticlePtr(Gluons(i),TmpGluons(counter))
!               counter=counter+1
!             enddo
!             rIn =1
!             rOut=0
!             do i=rIn,rOut
!               call CopyParticlePtr(Gluons(i),TmpGluons(counter))
!               counter=counter+1
!             enddo
!             tmp = cur_s_2s(TmpGluons(1:0),TmpScalar(1:1),(/0,0,0/) )
!             Res = Res + tmp
! 
! ! print *, "bg res",tmp
! ! print *, "ms res",Scalars(2)%Pol(1)*Scalars(3)%Pol(1)*Scalars(4)%Pol(1) * (0d0,1d0)/2d0/((Scalars(2)%Mom(1:Dv)+Scalars(3)%Mom(1:Dv)).ndot.(Scalars(2)%Mom(1:Dv)+Scalars(3)%Mom(1:Dv))) & 
! !                  * ((Scalars(2)%Mom(1:Dv)-Scalars(3)%Mom(1:Dv)).ndot.(Scalars(4)%Mom(1:Dv)+PMom1(1:Dv)))
! 




! NEW CODE
!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4).ne.0 ) print *, "wrong number of gluons in cur_s_4s"
!DEC$ ENDIF
      Res=(0d0,0d0)
      do n1a=0,NumGlu(1)
      do n3a=0,NumGlu(3)
      do n4a=0,NumGlu(4)
      do n1c=0,NumGlu(1)-n1a
      do n3c=0,NumGlu(3)-n3a
      do n4c=0,NumGlu(4)-n4a

         n1b = NumGlu(1)-n1a-n1c
         n3b = NumGlu(3)-n3a-n3c
         n4b = NumGlu(4)-n4a-n4c
         if( n1c.gt.0 .and. (n3c+n4c).gt.0  ) cycle
         if( n3c.gt.0 .and. (n1c+n4c).gt.0  ) cycle
         if( n4c.gt.0 .and. (n1c+n3c).gt.0  ) cycle

         rIn =n1a+n1c+1
         rOut=NumGlu(1)+NumGlu(2)+n3a
         Eps2 = cur_g_2s(Gluons(rIn:rOut),Scalars(2:3),(/1+n1b+NumGlu(2)+n3a,n1b,NumGlu(2),n3a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Scalars(2)%Mom + Scalars(3)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
         Eps2 = Eps2*PropFac1

            rIn =NumGlu(1)+NumGlu(2)+n3a+n3c+1
            rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
            Sca1 = cur_s_2s(Gluons(rIn:rOut),Scalars(4:4),(/n3b+n4a,n3b,n4a/) )! here was a bug
            PMom2(:) = Scalars(4)%Mom + SumMom(Gluons,rIn,rOut)! here was a bug
            if(n3b.ge.1 .or. n4a.ge.1) then
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Scalars(4)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Scalars(4)%Mass2).lt.PropCut ) cycle
               Sca1 = Sca1 * PropFac2
            endif

            if( n1c.gt.0 ) then
                rIn =n1a+1
                rOut=n1a+n1c
                EpsX(:) = cur_g(Gluons(rIn:rOut),1+n1c)
                if(n1c.gt.1) then
                    PMom4 = SumMom(Gluons,rIn,rOut)
                    PropFac4 = (0d0,-1d0)/sc_(PMom4,PMom4)
                    if( abs(sc_(PMom4,PMom4)).lt.PropCut ) cycle
                    EpsX(:) = EpsX(:)*PropFac4
                endif
                Sca0 = cggs(EpsX,Eps2) * Sca1
            elseif( n3c.gt.0) then
                rIn =NumGlu(1)+NumGlu(2)+n3a+1
                rOut=NumGlu(1)+NumGlu(2)+n3a+n3c
                EpsX(:) = cur_g(Gluons(rIn:rOut),1+n3c)
                if(n3c.gt.1) then
                    PMom4 = SumMom(Gluons,rIn,rOut)
                    PropFac4 = (0d0,-1d0)/sc_(PMom4,PMom4)
                    if( abs(sc_(PMom4,PMom4)).lt.PropCut ) cycle
                    EpsX(:) = EpsX(:)*PropFac4
                endif
                Sca0 = cggs(EpsX,Eps2) * Sca1
            elseif( n4c.gt.0 ) then
                rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+n4c
                EpsX(:) = cur_g(Gluons(rIn:rOut),1+n4c)
                if(n4c.gt.1) then
                    PMom4 = SumMom(Gluons,rIn,rOut)
                    PropFac4 = (0d0,-1d0)/sc_(PMom4,PMom4)
                    if( abs(sc_(PMom4,PMom4)).lt.PropCut ) cycle
                    EpsX(:) = EpsX(:)*PropFac4
                endif
                Sca0 = cgsg(Eps2,EpsX) * Sca1
            else
                if( Scalars(4)%PartType.gt.0 ) then
                    Sca0 = cgs(Eps2,PMom1,PMom1+PMom2) * Sca1
                else
                    Sca0 = cgbs(Eps2,PMom1,PMom1+PMom2) * Sca1
                endif
            endif

            PMom1 = Scalars(2)%Mom+Scalars(3)%Mom+Scalars(4)%Mom+SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+n4c)
            if(n1a.ge.1 .or. n4b.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Scalars(4)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Scalars(4)%Mass2).lt.PropCut ) cycle
               Sca0 = Sca0 * PropFac1
            endif


            TmpScalar(1)%Mom  => PMom1(:)
            TmpScalar(1)%Pol  => Sca0
            TmpScalar(1)%Mass => Scalars(4)%Mass
            TmpScalar(1)%Mass2=> Scalars(4)%Mass2
            TmpExtRef = -1
            TmpScalar(1)%ExtRef => TmpExtRef
            TmpScalar(1)%PartType => Scalars(4)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+n4c+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp = cur_s_2s(TmpGluons(1:counter-1),TmpScalar(1:1),(/counter-1,n1a,n4b/) )
            Res = Res + tmp

      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

RETURN

elseif( tag_f.eq.3 .and. NumGlu(2)+NumGlu(4).ne.0 ) then
    print *, "Error in cur_s_4s with tag_f=3; requested current is not available"
    stop
endif
endif




      Res=(0d0,0d0)
      do n1a=0,NumGlu(1)
      do n2a=0,NumGlu(2)
      do n4a=0,NumGlu(4)
      do n1c=0,NumGlu(1)-n1a
      do n2c=0,NumGlu(2)-n2a
      do n4c=0,NumGlu(4)-n4a

         n1b = NumGlu(1)-n1a-n1c
         n2b = NumGlu(2)-n2a-n2c
         n4b = NumGlu(4)-n4a-n4c
         if( n1c.gt.0 .and. (n2c+n4c).gt.0  ) cycle
         if( n2c.gt.0 .and. (n1c+n4c).gt.0  ) cycle
         if( n4c.gt.0 .and. (n1c+n2c).gt.0  ) cycle

         rIn =NumGlu(1)+n2a+n2c+1
         rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
         Eps2 = cur_g_2s(Gluons(rIn:rOut),Scalars(3:4),(/1+n2b+NumGlu(3)+n4a,n2b,NumGlu(3),n4a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Scalars(3)%Mom + Scalars(4)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
         Eps2 = Eps2*PropFac1

            rIn =n1a+n1c+1
            rOut=NumGlu(1)+n2a
            Sca1 = cur_s_2s(Gluons(rIn:rOut),Scalars(2:2),(/n2a+n1b,n1b,n2a/))
            PMom2(:) = Scalars(2)%Mom + SumMom(Gluons,rIn,rOut)
            if(n1b.ge.1 .or. n2a.ge.1) then
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Scalars(2)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Scalars(2)%Mass2).lt.PropCut ) cycle
               Sca1 = Sca1*PropFac2
            endif


            if( n1c.gt.0 ) then
                rIn =n1a+1
                rOut=n1a+n1c
                EpsX(:) = cur_g(Gluons(rIn:rOut),1+n1c)
                if(n1c.gt.1) then
                    PMom4 = SumMom(Gluons,rIn,rOut)
                    PropFac4 = (0d0,-1d0)/sc_(PMom4,PMom4)
                    if( abs(sc_(PMom4,PMom4)).lt.PropCut ) cycle
                    EpsX(:) = EpsX(:)*PropFac4
                endif
                Sca0 = cgsg(EpsX,Eps2) * Sca1!   here was a bug
            elseif( n2c.gt.0) then
                rIn =NumGlu(1)+n2a+1
                rOut=NumGlu(1)+n2a+n2c
                EpsX(:) = cur_g(Gluons(rIn:rOut),1+n2c)
                if(n2c.gt.1) then
                    PMom4 = SumMom(Gluons,rIn,rOut)
                    PropFac4 = (0d0,-1d0)/sc_(PMom4,PMom4)
                    if( abs(sc_(PMom4,PMom4)).lt.PropCut ) cycle
                    EpsX(:) = EpsX(:)*PropFac4
                endif
                Sca0 = csgg(Eps2,EpsX) * Sca1
            elseif( n4c.gt.0 ) then
                rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+n4c
                EpsX(:) = cur_g(Gluons(rIn:rOut),1+n4c)
                if(n4c.gt.1) then
                    PMom4 = SumMom(Gluons,rIn,rOut)
                    PropFac4 = (0d0,-1d0)/sc_(PMom4,PMom4)
                    if( abs(sc_(PMom4,PMom4)).lt.PropCut ) cycle
                    EpsX(:) = EpsX(:)*PropFac4
                endif
                Sca0 = csgg(Eps2,EpsX) * Sca1
            else
                  if( Scalars(2)%PartType.gt.0 ) then
                      Sca0 = csg(Eps2,PMom1,PMom2) * Sca1
                  else
                      Sca0 = cbsg(Eps2,PMom1,PMom2) * Sca1
                  endif
            endif


            PMom1 = Scalars(2)%Mom+Scalars(3)%Mom+Scalars(4)%Mom+SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+n4c)
            if(n1a.ge.1 .or. n4b.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Scalars(2)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Scalars(2)%Mass2).lt.PropCut ) cycle
               Sca0 = Sca0*PropFac1
            endif

            TmpScalar(1)%Mom  => PMom1(:)
            TmpScalar(1)%Pol  => Sca0
            TmpScalar(1)%Mass => Scalars(2)%Mass
            TmpScalar(1)%Mass2=> Scalars(2)%Mass2
            TmpExtRef = -1
            TmpScalar(1)%ExtRef => TmpExtRef
            TmpScalar(1)%PartType => Scalars(2)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+n4c+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp = cur_s_2s(TmpGluons(1:counter-1),TmpScalar(1:1),(/counter-1,n1a,n4b/) )
            Res = Res + tmp
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo



! adding the ssss vertex
      do n1a=0,NumGlu(1)
      do n2a=0,NumGlu(2)
      do n3a=0,NumGlu(3)
      do n4a=0,NumGlu(4)
         n1b = NumGlu(1)-n1a
         n2b = NumGlu(2)-n2a
         n3b = NumGlu(3)-n3a
         n4b = NumGlu(4)-n4a

            rIn =n1a+1
            rOut=NumGlu(1)+n2a
            Sca1 = cur_s_2s(Gluons(rIn:rOut),Scalars(2:2),(/n2a+n1b,n1b,n2a/))
            PMom1(:) = Scalars(2)%Mom + SumMom(Gluons,rIn,rOut)
            if(n1b.ge.1 .or. n2a.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Scalars(2)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Scalars(2)%Mass2).lt.PropCut ) cycle
               Sca1 = Sca1*PropFac1
            endif

            rIn =NumGlu(1)+n2a+1
            rOut=NumGlu(1)+NumGlu(2)+n3a
            Sca2 = cur_s_2s(Gluons(rIn:rOut),Scalars(3:3),(/n3a+n2b,n2b,n3a/))
            PMom2(:) = Scalars(3)%Mom + SumMom(Gluons,rIn,rOut)
            if(n2b.ge.1 .or. n3a.ge.1) then
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Scalars(3)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Scalars(3)%Mass2).lt.PropCut ) cycle
               Sca2 = Sca2*PropFac2
            endif

            rIn =NumGlu(1)+NumGlu(2)+n3a+1
            rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
            Sca3 = cur_s_2s(Gluons(rIn:rOut),Scalars(4:4),(/n4a+n3b,n3b,n4a/))
            PMom3(:) = Scalars(4)%Mom + SumMom(Gluons,rIn,rOut)
            if(n3b.ge.1 .or. n4a.ge.1) then
               PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3)-Scalars(4)%Mass2)
               if( abs(sc_(PMom3,PMom3)-Scalars(4)%Mass2).lt.PropCut ) cycle
               Sca3 = Sca3*PropFac3
            endif

            Sca0 = Sca1 * Sca2 * Sca3 * (0d0,-1d0)/2d0  

            PMom1 = Scalars(2)%Mom+Scalars(3)%Mom+Scalars(4)%Mom+SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a)
            if(n1a.ge.1 .or. n4b.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Scalars(2)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Scalars(2)%Mass2).lt.PropCut ) cycle
               Sca0 = Sca0*PropFac1
            endif

            TmpScalar(1)%Mom  => PMom1(:)
            TmpScalar(1)%Pol  => Sca0
            TmpScalar(1)%Mass => Scalars(2)%Mass
            TmpScalar(1)%Mass2=> Scalars(2)%Mass2
            TmpExtRef = -1
            TmpScalar(1)%ExtRef => TmpExtRef
            TmpScalar(1)%PartType => Scalars(2)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp = cur_s_2s(TmpGluons(1:counter-1),TmpScalar(1:1),(/counter-1,n1a,n4b/) )
            Res = Res + tmp
      enddo
      enddo
      enddo
      enddo



RETURN
END FUNCTION





!     assuming the color flow:   s = x--->----
      FUNCTION cur_s_2s_massCT(Gluons,Scalar,NumGlu) result(res)
      implicit none
      type(PtrToParticle) :: Gluons(1:),Scalar(2:2)      ! off-shell scalar is not included
      integer ::  NumGlu(0:2)
      integer :: ms1,m1,m2,ng1, ng2, ngluon
      complex(8) :: res,tmp 
      complex(8) :: sc2
      complex(8) :: k1(1:Dv),k2(1:Dv),k3(1:Dv)
      complex(8) :: e1(1:Dv)
      complex(8) :: e2(1:Dv)
      complex(8) :: k1sq,k2sq,k3sq



       ngluon = NumGlu(0)
       ng1 = NumGlu(1)           !#gluons to the left of a s-line
       ng2 = NumGlu(2)           !#gluons to the right of the s-line

       if (ngluon.ne.2) then
         call Error("scalar mass counterterm needs modification")
       endif
       res = (0d0,0d0)


if( NumGlu(1).eq.0 .and. NumGlu(2).eq.2 ) then

       k1 = Gluons(1)%Mom
       e1 = Gluons(1)%Pol

       sc2 = Scalar(2)%Pol(1)
       k2 = Scalar(2)%Mom
       if( Scalar(2)%PartType.gt.0 ) then
           tmp = csg(e1,k1,k2) * sc2
       else
           tmp = cbsg(e1,k1,k2) * sc2
       endif

       k2 = k1 + k2
       k2sq = sc_(k2,k2) - Scalar(2)%Mass2
       tmp =  (0d0,1d0)/k2sq*tmp! scalar prop
       ! here the CT is inserted
       tmp =  (0d0,1d0)/k2sq*tmp! scalar prop

       k1 = Gluons(2)%Mom
       e1 = Gluons(2)%Pol
       if( Scalar(2)%PartType.gt.0 ) then
           res = csg(e1,k1,k2) * tmp
       else
           res = cbsg(e1,k1,k2) * tmp
       endif


elseif( NumGlu(1).eq.2 .and. NumGlu(2).eq.0 ) then

       k1 = Gluons(2)%Mom
       e1 = Gluons(2)%Pol

       sc2 = Scalar(2)%Pol(1)
       k2 = Scalar(2)%Mom
       if( Scalar(2)%PartType.gt.0 ) then
           tmp = cgs(e1,k1,k2) * sc2
       else
           tmp = cgbs(e1,k1,k2) * sc2
       endif

       k2 = k1 + k2
       k2sq = sc_(k2,k2) - Scalar(2)%Mass2
       tmp =  (0d0,1d0)/k2sq*tmp! scalar prop
       ! here the CT is inserted
       tmp =  (0d0,1d0)/k2sq*tmp! scalar prop

       k1 = Gluons(1)%Mom
       e1 = Gluons(1)%Pol
       if( Scalar(2)%PartType.gt.0 ) then
           res = cgs(e1,k1,k2) * tmp
       else
           res = cgbs(e1,k1,k2) * tmp
       endif

elseif( NumGlu(1).eq.1 .and. NumGlu(2).eq.1 ) then

       k1 = Gluons(1)%Mom
       e1 = Gluons(1)%Pol

       sc2 = Scalar(2)%Pol(1)
       k2 = Scalar(2)%Mom
       if( Scalar(2)%PartType.gt.0 ) then
           tmp = cgs(e1,k1,k2) * sc2
       else
           tmp = cgbs(e1,k1,k2) * sc2
       endif

       k2 = k1 + k2
       k2sq = sc_(k2,k2) - Scalar(2)%Mass2
       tmp =  (0d0,1d0)/k2sq*tmp! scalar prop
       ! here the CT is inserted
       tmp =  (0d0,1d0)/k2sq*tmp! scalar prop

       k1 = Gluons(2)%Mom
       e1 = Gluons(2)%Pol
       if( Scalar(2)%PartType.gt.0 ) then
           res = csg(e1,k1,k2) * tmp
       else
           res = cbsg(e1,k1,k2) * tmp
       endif


!------

       k1 = Gluons(2)%Mom
       e1 = Gluons(2)%Pol

       sc2 = Scalar(2)%Pol(1)
       k2 = Scalar(2)%Mom
       if( Scalar(2)%PartType.gt.0 ) then
           tmp = csg(e1,k1,k2) * sc2
       else
           tmp = cbsg(e1,k1,k2) * sc2
       endif

       k2 = k1 + k2
       k2sq = sc_(k2,k2) - Scalar(2)%Mass2
       tmp =  (0d0,1d0)/k2sq*tmp! scalar prop
       ! here the CT is inserted
       tmp =  (0d0,1d0)/k2sq*tmp! scalar prop

       k1 = Gluons(1)%Mom
       e1 = Gluons(1)%Pol
       if( Scalar(2)%PartType.gt.0 ) then
           res = res + cgs(e1,k1,k2) * tmp
       else
           res = res + cgbs(e1,k1,k2) * tmp
       endif



endif



      end function cur_s_2s_massCT









END MODULE ModMyRecurrence


