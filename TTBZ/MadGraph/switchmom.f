C
C*********************************************************************************
C
C  Changes stuff for crossings
C
      subroutine switchmom( p1,p,ic,jc,nexternal )
      implicit none
      integer  nexternal
      integer  jc(nexternal), ic(nexternal)
      double precision  p1(0:3,nexternal), p(0:3,nexternal)
      integer  i,j
cc
      do i = 1,nexternal
         do j = 0,3
            p(j,ic(i)) = p1(j,i)
         enddo
      enddo
      do i = 1,nexternal
         jc(i) = 1
      enddo
      jc(ic(1)) = -1
      jc(ic(2)) = -1
ccc
      return
      end
