! write(*,*) "******************************************************"
! write(*,*) "Including Fact_Cuts"
! write(*,*) "******************************************************"
!DEC$ IF(_FactCheck .EQ.1)
   min = dmin1(dabs(MomExt(1:4,8).dot.MomExt(1:4,9)),& !g1-g2
               dabs(MomExt(1:4,6).dot.MomExt(1:4,8)),& !g1-q
               dabs(MomExt(1:4,7).dot.MomExt(1:4,8)),& !g1-qb
               dabs(MomExt(1:4,6).dot.MomExt(1:4,9)),& !1-4
               dabs(MomExt(1:4,7).dot.MomExt(1:4,9)),& !2-4
               dabs(MomExt(1:4,6).dot.MomExt(1:4,7)))/m_W**2 !1-2


  if(dabs(MomExt(1:4,8).dot.MomExt(1:4,9))/m_w**2.eq. min .and. min .lt.smin  ) then
     ! write(*,*) "1-recomb"
      if(dabs((MomExt(1:4,8)+MomExt(1:4,9)).dot.MomExt(1:4,6))/m_w**2.lt.smin ) then
          applyPSCut =.true.
           ! write(*,*) "1-1"
          goto 1276
      endif
      if(dabs((MomExt(1:4,8)+MomExt(1:4,9)).dot.MomExt(1:4,7))/m_w**2.lt.smin ) then
          applyPSCut =.true.
           ! write(*,*) "1-2"
          goto 1276
      endif
      if(dabs(MomExt(1:4,6).dot.MomExt(1:4,7))/m_w**2.lt.smin ) then
          applyPSCut =.true.
          ! write(*,*) "1-3"
          goto 1276
      endif
   endif

  if(dabs(MomExt(1:4,6).dot.MomExt(1:4,8))/m_w**2.eq. min .and. min .lt.smin  ) then
     ! write(*,*) "2-recomb"
      if(dabs((MomExt(1:4,6)+MomExt(1:4,8)).dot.MomExt(1:4,9))/m_w**2.lt.smin ) then
          applyPSCut =.true.
           ! write(*,*) "2-1"
          goto 1276
      endif
      if(dabs((MomExt(1:4,6)+MomExt(1:4,8)).dot.MomExt(1:4,7))/m_w**2.lt.smin ) then
          applyPSCut =.true.
           ! write(*,*) "2-2"
          goto 1276
      endif
      if(dabs(MomExt(1:4,9).dot.MomExt(1:4,7))/m_w**2.lt.smin ) then
           ! write(*,*) "2-3"
          applyPSCut =.true.
          goto 1276
      endif
  endif
  if(dabs(MomExt(1:4,6).dot.MomExt(1:4,9))/m_w**2.eq. min .and. min .lt.smin  ) then
     ! write(*,*) "3-recomb"
      if(dabs((MomExt(1:4,6)+MomExt(1:4,9)).dot.MomExt(1:4,8))/m_w**2.lt.smin ) then
           ! write(*,*) "3-1"
          applyPSCut =.true.
          goto 1276
      endif
      if(dabs((MomExt(1:4,6)+MomExt(1:4,9)).dot.MomExt(1:4,7))/m_w**2.lt.smin ) then
           ! write(*,*) "3-2"
          applyPSCut =.true.
          goto 1276
      endif
      if(dabs(MomExt(1:4,8).dot.MomExt(1:4,7))/m_w**2.lt.smin ) then
           ! write(*,*) "3-3"
          applyPSCut =.true.
          goto 1276
      endif
  endif
  if(dabs(MomExt(1:4,7).dot.MomExt(1:4,8))/m_w**2.eq. min .and. min .lt.smin  ) then
     ! write(*,*) "4-recomb"
      if(dabs((MomExt(1:4,7)+MomExt(1:4,8)).dot.MomExt(1:4,9))/m_w**2.lt.smin ) then
           ! write(*,*) "4-1"
          applyPSCut =.true.
          goto 1276
      endif
      if(dabs((MomExt(1:4,7)+MomExt(1:4,8)).dot.MomExt(1:4,6))/m_w**2.lt.smin ) then
           ! write(*,*) "4-2"
          applyPSCut =.true.
          goto 1276
      endif
      if(dabs(MomExt(1:4,9).dot.MomExt(1:4,6))/m_w**2.lt.smin ) then
           ! write(*,*) "4-3"
          applyPSCut =.true.
          goto 1276
      endif
  endif
  if(dabs(MomExt(1:4,7).dot.MomExt(1:4,9))/m_w**2.eq. min .and. min .lt.smin   ) then
     ! write(*,*) "5-recomb"
      if(dabs((MomExt(1:4,7)+MomExt(1:4,9)).dot.MomExt(1:4,8))/m_w**2.lt.smin ) then
           ! write(*,*) "5-1"
          applyPSCut =.true.
          goto 1276
      endif
      if(dabs((MomExt(1:4,7)+MomExt(1:4,9)).dot.MomExt(1:4,6))/m_w**2.lt.smin ) then
           ! write(*,*) "5-2"
          applyPSCut =.true.
          goto 1276
      endif
      if(dabs(MomExt(1:4,8).dot.MomExt(1:4,6))/m_w**2.lt.smin ) then
           ! write(*,*) "5-3"
          applyPSCut =.true.
           goto 1276
      endif
  endif
  if(dabs(MomExt(1:4,6).dot.MomExt(1:4,7))/m_w**2.eq. min .and. min .lt.smin ) then
     ! write(*,*) "6-recomb"
      if(dabs((MomExt(1:4,7)+MomExt(1:4,6)).dot.MomExt(1:4,8))/m_w**2.lt.smin ) then
           ! write(*,*) "6-1"
          applyPSCut =.true.
          goto 1276
      endif
      if(dabs((MomExt(1:4,7)+MomExt(1:4,6)).dot.MomExt(1:4,9))/m_w**2.lt.smin ) then
           ! write(*,*) "6-2"
          applyPSCut =.true.
          goto 1276
      endif
      if(dabs(MomExt(1:4,8).dot.MomExt(1:4,9))/m_w**2.lt.smin ) then
           ! write(*,*) "6-3"
          applyPSCut =.true.
           goto 1276
      endif
  endif

! no recombination
!  If ((dabs(MomExt(1:4,9).dot.MomExt(1:4,8))/m_W**2 .lt. smin ) .or. (dabs(MomExt(1:4,7).dot.MomExt(1:4,9))/m_W**2 .lt. smin )) then
!!           ! write(*,*) "6"
!     applyPSCut =.true.
!     goto 1276
!  endif
1276 continue
!DEC$ ENDIF

