! Module IntegerPartition
! PartitionFunction(n):         returns number of inequivalent partitions of the integer 'n'
! GetIntPartition(n,Partition): returns array 'Partition(:,:)' which contains all partitions of the integer 'n'
!                               'Partition(NPart,Col)' has to be allocated before the call, 'NPart'=number of partitions, 'Col'=elements of a partition
!
! routines are not optimized for speed


MODULE ModIntegerPartition
implicit none
save

public  :: PartitionFunction, GetIntPartition
private :: AppendToList, CalcIntPartition

CONTAINS




RECURSIVE FUNCTION PartitionFunction(n) result(res)
implicit none
integer :: n
integer :: k,res

   if( n.eq.0 ) then
      res = 1
      return
   endif

   res = 0
   do k=1,n
      res = res + (-1)**(k+1) * ( PartitionFunction(n-(3*k*k-k)/2) + PartitionFunction(n-(3*k*k+k)/2) )
   enddo
return
END FUNCTION



RECURSIVE SUBROUTINE CalcIntPartition(n,limit,TheList,c,Partition)
implicit none
integer :: n,limit,c
integer, allocatable :: TheList(:)
integer, allocatable :: TmpList(:)
integer, allocatable :: Partition(:,:)
integer :: i

   allocate(TmpList(1:size(TheList)))
   if( n.gt.0 ) then
      do i=min(n,limit),1,-1
         TmpList(:) = AppendToList(TheList(:),i)
         call CalcIntPartition(n-i,i,TmpList,c,Partition)
      enddo
   else
      Partition(c,:) = TheList(:)
      c=c+1
   endif
return
END SUBROUTINE



FUNCTION AppendToList(Partition,Element) 
implicit none
integer :: Partition(:)
integer :: AppendToList(1:size(Partition))
integer :: i,Element

   do i=1,size(Partition)
      AppendToList(i) = Partition(i)
      if(Partition(i).eq.0) then
         AppendToList(i) = Element
	 AppendToList(i+1:size(Partition))=0
	 return
      endif
   enddo
return
END FUNCTION

 

SUBROUTINE GetIntPartition(n,Partition)
implicit none
integer :: n
integer, allocatable :: Partition(:,:)
integer :: j,i,Counter,AllocStatus
integer,allocatable :: DummyList(:)


   do i=1,size(Partition,dim=1)
      do j=1,size(Partition,dim=2)
         Partition(i,j) = 0
      enddo
   enddo
   
   if( n.le.0 ) return
   
   allocate(DummyList(1:n),stat=AllocStatus)
   if( AllocStatus .ne. 0 ) print *,"Error in memory allocation for DummyList(1:n)"
   do j=1,n
      DummyList(j) = 0
   enddo

   Counter=1
   call CalcIntPartition(n,n,DummyList,Counter,Partition)

return
END SUBROUTINE


END MODULE
