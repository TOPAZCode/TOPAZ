! Module Permutations
! GetPermutation(Set,PermSet): returns array 'PermSet(:,:)' which contains all permutations of the set 'Set(:)'
!                             'PermSet(NPerm,Col)' has to be allocated before the call, 'NPerm'=number of permutations, 'Col'=elements of the permutation
!
! routines are not optimized for speed


MODULE ModPermutations
implicit none
save

public  :: GetPermutation
private :: CalcPermutation,Remove,Factorial

CONTAINS




RECURSIVE SUBROUTINE CalcPermutation(TheSet,Perm,PermSet,L,N,c,NoDblCount)
implicit none
integer, allocatable :: TheSet(:),Perm(:),TheList(:),PermSet(:,:)
logical, OPTIONAL :: NoDblCount
integer :: L,N,i,AllocStatus,c

   if(L.le.N) then
      allocate(TheList(size(TheSet)-1),stat=AllocStatus)   
      if( AllocStatus .ne. 0 ) print *,"Error in memory allocation for TheList(1:n)"
      do i=1,size(TheSet)
         Perm(L) = TheSet(i);
	 TheList(1:size(TheSet)-1) = Remove(TheSet,i)
	 call CalcPermutation(TheList,Perm,PermSet,L+1,N,c,NoDblCount)
      enddo
   else
      if( present(NoDblCount) ) then
         if( NoDblCount ) then
	    call InsertIfNew(PermSet,Perm,c)
	 else
	    PermSet(c,1:N) = Perm(1:N)
            c=c+1
	 endif
      else
         PermSet(c,1:N) = Perm(1:N)
         c=c+1
      endif
   endif
return
END SUBROUTINE



SUBROUTINE InsertIfNew(PermSet,Perm,c)
implicit none
integer, allocatable :: Perm(:),PermSet(:,:)
integer :: c,n,per,col
logical :: EqualPerm

   n = size(Perm)
   do per=1,c-1
      EqualPerm=.true.
      do col=1,n
         if( PermSet(per,col).ne.Perm(col) ) EqualPerm=.false.
      enddo
      if(EqualPerm) return
   enddo

   PermSet(c,1:N) = Perm(1:N)
   c=c+1
return
END SUBROUTINE



FUNCTION Remove(List,Element) 
implicit none
integer :: List(:),Element,j1,j2
integer :: Remove(1:size(List)-1)

      j1=1; j2=1;      
      do while( j2.le.size(List) )
	    if( j2.eq.Element ) then
	    	j2=j2+1
		cycle
            endif
	    Remove(j1) = List(j2)
	    j1=j1+1
	    j2=j2+1
      enddo
return
END FUNCTION



FUNCTION Factorial(N)
implicit none
integer Factorial,n,i

   factorial = 1
   do i=n,1,-1
       factorial = factorial * i
   end do
return
END FUNCTION




SUBROUTINE GetPermutation(Set,PermSet,NoDblCount)
implicit none
integer, allocatable :: Set(:),PermSet(:,:),PermVec(:)
logical, OPTIONAL :: NoDblCount
integer N,AllocStatus,c

   N=size(Set)
   allocate(PermVec(1:N),stat=AllocStatus)
   if( AllocStatus .ne. 0 ) print *,"Error in memory allocation for PermVec(1:n)"
   
   c=1
   if( present(NoDblCount) ) then 
      call CalcPermutation(Set,PermVec,PermSet,1,N,c,NoDblCount)
   else
      call CalcPermutation(Set,PermVec,PermSet,1,N,c)
   endif
return
END SUBROUTINE


END MODULE
