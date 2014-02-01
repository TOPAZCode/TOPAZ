! Compile:  ifort -o BinnedLogL BinnedLogL.f90
! Compile:  mpif90 -D_UseMPI=1 -o BinnedLogL BinnedLogL.f90
! Compile:  safety options:  -O0 -fpp -implicitnone -zero -check bounds -check pointer -warn interfaces -ftrapuv
program BinnedLogL
  use ifport
  implicit none
  integer :: NumArgs,Histo,NPseudoExp,NPseudoExp_Worker
  real(8) :: PreFactor,Data
  character :: operation*(10),H0_infile*(50),H1_infile*(50),Data_str*(5),Histo_str*(5),NPseudoExp_str*(9),dummy*(1),outfile*(50),DeltaN_str*(5),filename_out1*(50),filename_out2*(50)
  character(len=*),parameter :: fmt1 = "(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)"
character(len=*),parameter :: fmt2 = "(I5,A,2X,1PE14.7,A,2X,I9,A,2X,1PE14.7,A,2X,I9,A,2X,1PE14.7,A,2X,1PE14.7,A,2X,1PE14.7,A)"
  integer,parameter :: MaxBins=1100, MaxEvents=100000,MaxExp=1000000
  integer,parameter :: NumLLbins=1100
  real(8) :: BinVal(1:2,1:MaxBins)=-1d-99,Value(1:2,1:MaxBins)=-1d-99,Error(1:2,1:MaxBins)=-1d-99,BinSize(1:2)
  integer :: NHisto(1:2)=-999999,Hits(1:2,1:MaxBins)=-999999,NumBins0=0,NumBins1=0,NumBins,iPseudoExp,iBin,iHypothesis,LLBin,NumExpectedEvents(1:2)
  real(8) :: PoissonMax(1:2,1:MaxBins),LLRatio
  logical :: GotNumEvents(1:MaxBins)
  integer :: ExpectedEvents(1:2,1:MaxBins),TryEvts,ObsEvents(1:2,1:MaxBins),PlotObsEvts(1:2,1:10000)
  real(8) :: nran(1:2),offset,sran,DeltaN
  real(8) :: sigmatot(1:2),check(1:2),alpha(1:2,1:MaxBins),IntLLRatio(1:2),checkalpha(1:MaxBins)
  real(8) :: alphamin,alphamax,betamin,betamax,rescale(1:2),alpha_ave,beta_ave,sigs
  real(8) :: LLRatio_array(1:2,1:MaxExp),LLRatio_array_tmp(1:2,1:MaxExp),LLRatio_min,LLRatio_max,WhichBin
  integer :: j,i,s,SUA,MPI_Rank,worker,MPI_NUM_PROCS,ierror
  type :: Histogram
     integer :: NBins
     real(8) :: BinSize
     real(8) :: LowVal
     integer :: Hits(1:1000)
  end type Histogram
  type(Histogram) :: LLHisto(1:2)
!DEC$ IF(_UseMPI .EQ.1)
  include 'mpif.h'
  integer status(MPI_STATUS_SIZE)
     call MPI_INIT(ierror)
     call MPI_COMM_RANK(MPI_COMM_WORLD,MPI_Rank,ierror)
     call MPI_COMM_SIZE (MPI_COMM_WORLD,MPI_NUM_PROCS,ierror)
!DEC$ ELSE
     MPI_Rank=0
     MPI_NUM_PROCS=1
!DEC$ ENDIF



  LLHisto(1)%NBins   = 1000
  LLHisto(1)%BinSize = 1d0
  LLHisto(1)%LowVal  = -500d0
  LLHisto(1)%Hits(:) = 0
  
  LLHisto(2)%NBins   =  LLHisto(1)%NBins   
  LLHisto(2)%BinSize =  LLHisto(1)%BinSize 
  LLHisto(2)%LowVal  =  LLHisto(1)%LowVal  
  LLHisto(2)%Hits(:) =  LLHisto(1)%Hits(:) 

  NumArgs=NArgs()-1
  PreFactor=8d0
  SUA=2


  if (SUA .ne. 1 .and. SUA .ne. 2 .and. DeltaN .ne. 1d0) then
     print *, "WARNING : SCALE UNCERTAINTY NOT USED!"
     stop
  endif

  call GetArg(1,operation)
  
  if (trim(operation) .ne. 'GetLogL') then
     print *, 'only available operation is GetLogL'
     stop
  endif

  call GetArg(2,H0_infile)
  call GetArg(3,H1_infile)
  call GetArg(4,Histo_str)
  call GetArg(5,Data_str)
  call GetArg(6,NPseudoExp_str)
  call GetArg(7,DeltaN_str)
  call GetArg(8,outfile)
  filename_out1=trim(outfile)//".dat"
  filename_out2=trim(outfile)//".out"
  read(Histo_str,"(I2)") Histo
  read(Data_str,"(F10.6)") Data
  read(NPseudoExp_str,"(I9)") NPseudoExp
  read(DeltaN_str,"(F9.6)") DeltaN
  write(*,*) "Performing binned log-likelihood analysis."
  write(*,"(A,I2,A,F10.6,A,I9,A,F9.6)") "Reading histogram ",Histo," for likelihood analysis with Lumi=",Data,"fb^-1 using ",NPseudoExp," pseudo-experiments, and scale uncertainty", DeltaN
  write(*,*) 'Using input files:', trim(H0_infile),trim(H1_infile)
  write(*,*) "Writing output to files : ", trim(filename_out1), " and ", trim(filename_out2)
  if (NPseudoExp .gt. MaxExp) then
     print *, "Too many experiments!"
     print *, MaxExp, NPseudoExp
     stop
  endif
  sigmatot(1:2) = 0d0
  NumBins0=0
  NumBins1=0
!--------------------------------------------------
!          1. reading input files
!--------------------------------------------------
  open(unit=112,file=trim(H0_infile),form='formatted',access='sequential')  ! open input file1
  open(unit=113,file=trim(H1_infile),form='formatted',access='sequential')  ! open input file
  open(unit=12,file=trim(H0_infile),form='formatted',access='sequential')  ! open input file1
  open(unit=13,file=trim(H1_infile),form='formatted',access='sequential')  ! open input file
  open(unit=14,file=trim(filename_out1),form='formatted',access='sequential')  ! open output file
  open(unit=15,file=trim(filename_out2),form='formatted',access='sequential')  ! open output file

  write(14,*) "Performing binned log-likelihood analysis."
  write(14,"(A,I2,A,F10.6,A,I9,A,F9.6)") "Reading histogram ",Histo," for likelihood analysis with Lumi=",Data,"fb^-1 using ",NPseudoExp," pseudo-experiments, and scale uncertainty", DeltaN
  write(14,*) 'Using input files:', trim(H0_infile)," and  ", trim(H1_infile)
  write(14,*) "Writing output to files : ", trim(filename_out1), " and ", trim(filename_out2)


!-------------------------------------------------------------------------
! 1a. Because our Delta phi_ll histogram excludes a small bin at the end,
! first use histo #1 to get the total cross-section
!-------------------------------------------------------------------------

  do while(.not.eof(112))!   reading input file 1
     read(unit=112,fmt="(A)") dummy
     if(dummy(1:1).eq."#") cycle
     backspace(unit=112) ! go to the beginning of the line
     
     read(unit=112,fmt=fmt1) NHisto(1),dummy,BinVal(1,NumBins0+1),dummy,Value(1,NumBins0+1),dummy,Error(1,NumBins0+1),dummy,Hits(1,NumBins0+1),dummy

!    if( NHisto(1).ne.1 ) cycle
    if( NHisto(1).ne.Histo ) cycle
     NumBins0=NumBins0 + 1
  enddo
  
  do while(.not.eof(113))!   reading input file 2
     read(unit=113,fmt="(A)") dummy
     if(dummy(1:1).eq."#") cycle
     backspace(unit=113) ! go to the beginning of the line
     
     read(unit=113,fmt=fmt1) NHisto(2),dummy,BinVal(2,NumBins1+1),dummy,Value(2,NumBins1+1),dummy,Error(2,NumBins1+1),dummy,Hits(2,NumBins1+1),dummy
!    if( NHisto(2).ne.1 ) cycle
    if( NHisto(2).ne.Histo ) cycle
     NumBins1=NumBins1 + 1
  enddo

  BinSize(1) = BinVal(1,2)-BinVal(1,1)
  BinSize(2) = BinVal(2,2)-BinVal(2,1)
  do iBin=1,NumBins1! calculate total cross section
       sigmatot(1) = sigmatot(1) + Value(1,iBin) * BinSize(1)
       sigmatot(2) = sigmatot(2) + Value(2,iBin) * BinSize(2)
  enddo

! now reset everything to zero
  NumBins0=0
  NumBins1=0
  BinSize=0d0
  NHisto=0
  BinVal=0d0
  Value=0d0
  Error=0d0
  Hits=0

!------------------------------------------------------------------------- 
! 1b. Now read in the actual distribution of interest
!-------------------------------------------------------------------------

  do while(.not.eof(12))!   reading input file 1
     read(unit=12,fmt="(A)") dummy
     if(dummy(1:1).eq."#") cycle
     backspace(unit=12) ! go to the beginning of the line
     read(unit=12,fmt=fmt1) NHisto(1),dummy,BinVal(1,NumBins0+1),dummy,Value(1,NumBins0+1),dummy,Error(1,NumBins0+1),dummy,Hits(1,NumBins0+1),dummy
    if( NHisto(1).ne.Histo ) cycle
     NumBins0=NumBins0 + 1
  enddo
  
  do while(.not.eof(13))!   reading input file 2
     read(unit=13,fmt="(A)") dummy
     if(dummy(1:1).eq."#") cycle
     backspace(unit=13) ! go to the beginning of the line
     
     read(unit=13,fmt=fmt1) NHisto(2),dummy,BinVal(2,NumBins1+1),dummy,Value(2,NumBins1+1),dummy,Error(2,NumBins1+1),dummy,Hits(2,NumBins1+1),dummy
     if( NHisto(2).ne.Histo ) cycle
     NumBins1=NumBins1 + 1
  enddo


  if( NumBins0 .ne. NumBins1 ) then
     print *, "Error: Number of bins in input file 1 and 2 are different: ",NumBins0,NumBins1
     stop
  else
     NumBins=NumBins0
  endif
  
  BinSize(1) = BinVal(1,2)-BinVal(1,1)
  BinSize(2) = BinVal(2,2)-BinVal(2,1)


  NumExpectedEvents(1) = int(Data * sigmatot(1) * PreFactor) ! factor of 8 from lepton species
  NumExpectedEvents(2) = int(Data * sigmatot(2) * PreFactor) 
  write(*,"(A,1PE16.8,A,I6)") "Total cross section of input file 1: ",sigmatot(1),"   <-->   Number of events: ",NumExpectedEvents(1)
  write(*,"(A,1PE16.8,A,I6)") "Total cross section of input file 2: ",sigmatot(2),"   <-->   Number of events: ",NumExpectedEvents(2)
  write(14,"(A,1PE16.8,A,I6)") "Total cross section of input file 1: ",sigmatot(1),"   <-->   Number of events: ",NumExpectedEvents(1)
  write(14,"(A,1PE16.8,A,I6)") "Total cross section of input file 2: ",sigmatot(2),"   <-->   Number of events: ",NumExpectedEvents(2)

!  check(1:2) = 0d0
!  do iBin=1,NumBins! check
!      check(1) = check(1) + Value(1,iBin)
!      check(2) = check(2) + Value(2,iBin)
!  enddo
!  write(*,"(A,1PE16.8)") "Sum of bins in input file 1: ",check(1)
!  write(*,"(A,1PE16.8)") "Sum of bins in input file 2: ",check(2)
!  write(*,*) ""

! print the input histograms 
  write(*,"(2X,A,16X,A,11X,A,16X,A)") "NBin|","Input file 1","|","Input file 2"
  do iBin=1,NumBins
     write(*,fmt="(2X,1I3,A,2X,1PE10.3,2X,1PE23.16,A,2X,1PE10.3,2X,1PE23.16)") iBin," | ",BinVal(1,iBin),Value(1,iBin)," | ",BinVal(2,iBin),Value(2,iBin)
     if( dabs(BinVal(1,iBin)-BinVal(2,iBin)).gt.1d-6 ) then
        print *, "Error: Different bin sizes in input files 1 and 2"
        stop
     endif
  enddo
  write(14,*) ""
  
! -----------------------------------------------------------------
! 2. Find the expected values for null and alt. hypothesis in each bin
! -----------------------------------------------------------------

! Null hypothesis first

  do iHypothesis=1,2    ! 1=null, 2=alt
     do iBin=1,NumBins
!     ExpectedEvents0(iBin)=int(Value(1,iBin)*Data*PreFactor*BinSize0)
!     ExpectedEvents1(iBin)=int(Value(2,iBin)*Data*PreFactor*BinSize1)
!
!     ! now find the value of the Poisson distribution at this maximum
!     PoissonMax0(iBin)=Poisson(ExpectedEvents0(iBin),ExpectedEvents0(iBin))
!     PoissonMax1(iBin)=Poisson(ExpectedEvents1(iBin),ExpectedEvents1(iBin))
        ExpectedEvents(iHypothesis,iBin)=int(Value(iHypothesis,iBin)*Data*PreFactor*BinSize(iHypothesis))
        if (ExpectedEvents(iHypothesis,iBin).eq.0) then
           PoissonMax(iHypothesis,iBin)=0d0
        else
           PoissonMax(iHypothesis,iBin)=Poisson(ExpectedEvents(iHypothesis,iBin),ExpectedEvents(iHypothesis,iBin))
        endif
     enddo
  enddo
  
!  call random_seed()
  call init_random_seed()
! -----------------------------------------------------------------
! 3.1 Generate Poisson distribution about expected null value in each bin
! -----------------------------------------------------------------
  NPseudoExp_Worker = NPseudoExp/(MPI_NUM_PROCS)
  print *, "MPI rank ",MPI_Rank," running ",NPseudoExp_Worker," pseudoexperiments"
!DEC$ IF(_UseMPI .EQ.1)
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
!DEC$ ENDIF
!  call random_seed()
  PlotObsEvts=0
  do iHypothesis=1,2
!           call random_seed()
     print *, 'Hypothesis ', iHypothesis
     do iPseudoExp=1,NPseudoExp_Worker
!        if( mod(iPseudoExp,10000).eq.0 ) print *, "Pseudo experiment ",iPseudoExp,"/",NPseudoExp
        GotNumEvents=.false.
        do iBin=1,NumBins
           do while (.not. GotNumEvents(iBin))
              call random_number(nran(1:2))
              nran(1)=nran(1)*PoissonMax(iHypothesis,iBin)        
              TryEvts=int(5d0*ExpectedEvents(iHypothesis,iBin)*nran(2))
              if (Poisson(ExpectedEvents(iHypothesis,iBin),TryEvts) .gt. nran(1)) then
                 ObsEvents(iHypothesis,iBin)=TryEvts
                 GotNumEvents(iBin)=.true.
              endif
              if (ExpectedEvents(iHypothesis,iBin) .eq. 0) then
                 ObsEvents(iHypothesis,iBin)=0
                 GotNumEvents(iBin)=.true.
              endif
           enddo
           if (SUA .eq. 1) then
              call random_number(sran)
              sran=sran*(DeltaN**2-1d0)/DeltaN + 1d0/DeltaN
              ObsEvents(iHypothesis,iBin)=ObsEvents(iHypothesis,iBin)*sran
           endif

           if (iBin .eq. 1 .and. ObsEvents(iHypothesis,iBin) .gt. 0) then
              PlotObsEvts(iHypothesis,ObsEvents(iHypothesis,iBin))=PlotObsEvts(iHypothesis,ObsEvents(iHypothesis,iBin))+1
           endif

        enddo

! -----------------------------------------------------------------! 3.1a  change number of observed events in all bins for scale uncertainty    
! -----------------------------------------------------------------
        if (SUA .eq. 2) then
           call random_number(sran)
           sran=sran*(DeltaN**2-1d0)/DeltaN + 1d0/DeltaN
           write(201,*) sran
           ObsEvents(iHypothesis,1:NumBins)=ObsEvents(iHypothesis,1:NumBins)*sran
        endif
! -----------------------------------------------------------------
! 3.2 Now find the log likelihood, with a Poisson distr in each bin
! -----------------------------------------------------------------

     LLRatio=0d0
     do iBin=1,NumBins
        if (ObsEvents(iHypothesis,iBin) .ne. 0 .and. ExpectedEvents(1,iBin) .ne. 0 .and. ExpectedEvents(2,iBin) .ne. 0) then
           LLRatio=LLRatio+ObsEvents(iHypothesis,iBin)*dlog(1d0*ExpectedEvents(1,iBin)/ExpectedEvents(2,iBin))
!           if (iPseudoExp .eq. 1) then
           write(301,*) iBin, ObsEvents(iHypothesis,iBin),dlog(1d0*ExpectedEvents(1,iBin)/ExpectedEvents(2,iBin)),ObsEvents(iHypothesis,iBin)*dlog(1d0*ExpectedEvents(1,iBin)/ExpectedEvents(2,iBin)),LLRatio
!        endif
        endif
     enddo

! offset to get the distributions with the histogram limits
     if (iHypothesis .eq. 1 .and. iPseudoExp .eq. 1) then
        offset=-100d0*(int(LLRatio)/100)
        print *, 'using offset of ', offset
     endif
     
     LLRatio=LLRatio+offset

     if (iPseudoExp .eq. 1) then
        print *, 'LLRatio=',iHypothesis,LLRatio
     endif
     LLRatio_array(iHypothesis,iPseudoExp)=LLRatio
        
  enddo!  iPseudoExp

  do i=1,1000
     s=101+iHypothesis
     write(s,*) i,PlotObsEvts(iHypothesis,i)
  enddo

enddo!   iHypothesis

  print *, "MPI rank ",MPI_Rank," finished"
!DEC$ IF(_UseMPI .EQ.1)
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
   if( MPI_Rank.eq.0 ) then
      do worker=1,MPI_NUM_PROCS-1
        call MPI_Recv(LLRatio_array_tmp, 2*MaxExp,MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE ,1,MPI_COMM_WORLD, status,ierror)
        print *, "received array from worker",worker
        LLRatio_array(1,NPseudoExp_Worker*worker+1:NPseudoExp_Worker*(worker+1)) = LLRatio_array_tmp(1,1:NPseudoExp_Worker)
        LLRatio_array(2,NPseudoExp_Worker*worker+1:NPseudoExp_Worker*(worker+1)) = LLRatio_array_tmp(2,1:NPseudoExp_Worker)
      enddo
   elseif( MPI_Rank.gt.0 ) then
        call MPI_Send(LLRatio_array, 2*MaxExp,MPI_DOUBLE_PRECISION, 0 ,1,MPI_COMM_WORLD,ierror)
   endif
!DEC$ ENDIF


if( MPI_Rank.eq.0 ) then 

!     if (iHypothesis .eq. 1 .and. iPseudoExp .eq. 1) then
!        print *, LLRatio,offset
!
!        if (LLRatio .lt. LLHisto1%LowVal .or. LLRatio .gt. (LLHisto1%LowVal+LLHisto1%NBins*LLHisto1%BinSize)) then
!            print *, 'LLRatio out  of binning range! Reset binning and rerun'
!            print *, LLRatio
!            stop
!         endif
!      endif

! -----------------------------------------------------------------
! 3.3 Bin the log likelihood
! -----------------------------------------------------------------

!      WhichBin = (LLRatio-LLHisto1%LowVal)/LLHisto1%BinSize + 1
!      if( WhichBin.lt.0 ) WhichBin = 1
!      if( WhichBin.gt.LLHisto1%NBins ) WhichBin = LLHisto1%NBins
!      LLHisto1%Hits(WhichBin) = LLHisto1%Hits(WhichBin) + 1
      

 !  enddo! iPseudoExp


!   do iPseudoExp=1,1000    ! this is just a check 
!      write(201,*) iPseudoExp,PlotObsEvts(1,iPseudoExp)
!   enddo
   
   LLRatio_max=-1d-6
   LLRatio_min=1d6
   
   do iHypothesis=1,2
      do iPseudoExp=1,NPseudoExp
         if (LLRatio_array(iHypothesis,iPseudoExp) .gt. LLRatio_max) then
            LLRatio_max=LLRatio_array(iHypothesis,iPseudoExp) 
         endif
         if (LLRatio_array(iHypothesis,iPseudoExp) .lt. LLRatio_min) then
            LLRatio_min=LLRatio_array(iHypothesis,iPseudoExp) 
         endif
      enddo
   enddo

   LLRatio_max=ceiling(LLRatio_max)
   LLRatio_min=floor(LLRatio_min)
   
   if (LLRatio_max .lt. LLRatio_min) then
      print *, 'ERROR: MAX OF LL RATIO IS SMALLER THAN MIN!'
      print *, LLRatio_max,LLRatio_min
      stop
   endif
   if (LLRatio_min .eq. LLRatio_max) then
      LLRatio_max=LLRatio_max+1d0
      LLRatio_min=LLRatio_min-1d0
   endif
   print *, LLRatio_min,LLRatio_max
   
   LLHisto(1)%NBins   = 1000
   LLHisto(1)%BinSize = (LLRatio_max-LLRatio_min)/LLHisto(1)%NBins
   LLHisto(1)%LowVal  = LLRatio_min
   LLHisto(1)%Hits(:) = 0
   print *, 'using ', LLHisto(1)%BinSize, ' bins in range ',LLRatio_min, LLRatio_max
   LLHisto(2)%NBins   = LLHisto(1)%NBins
   LLHisto(2)%BinSize = LLHisto(1)%BinSize
   LLHisto(2)%LowVal  = LLHisto(1)%LowVal
   LLHisto(2)%Hits(:) = LLHisto(1)%Hits(:)
   
   do iHypothesis=1,2
      do iPseudoExp=1,NPseudoExp
         LLRatio=LLRatio_array(iHypothesis,iPseudoExp)
         WhichBin = (LLRatio-LLHisto(iHypothesis)%LowVal)/LLHisto(iHypothesis)%BinSize + 1
         WhichBin=int(WhichBin)
         if( WhichBin.lt.0 ) WhichBin = 1
         if( WhichBin.gt.LLHisto(iHypothesis)%NBins ) WhichBin = LLHisto(iHypothesis)%NBins
         LLHisto(iHypothesis)%Hits(WhichBin) = LLHisto(iHypothesis)%Hits(WhichBin) + 1
      enddo

! -----------------------------------------------------------------
! 4. Calculate the integral of the log likelihood ratio curve, at each bin
! -----------------------------------------------------------------


      do LLbin=1,LLHisto(iHypothesis)%NBins
         IntLLRatio(iHypothesis)=IntLLRatio(iHypothesis)+&
              LLHisto(iHypothesis)%BinSize * LLHisto(iHypothesis)%Hits(LLbin)
         alpha(iHypothesis,LLBin)=IntLLRatio(iHypothesis)/(NPseudoExp*LLHisto(iHypothesis)%BinSize)
      enddo
   enddo

!  do LLbin=1,LLHisto1%NBins
!     IntH0=IntH0+LLHisto1%BinSize*LLHisto1%Hits(LLbin)
!     alpha(LLbin)=IntH0/(NPseudoExp*LLHisto1%BinSize)
!     print *, LLbin, LLHisto1%LowVal + LLbin*LLHisto1%BinSize, LLHisto1%Hits(LLbin), alpha(LLBin)
!  enddo

   print *, ""
   print *, "Writing log-likelihood distribution to output file"
   print *, ""
   do LLbin=1,LLHisto(1)%NBins
      !         print *, LLBin
      write(14,"(2X,I4,2X,1PE16.8,2X,I10,2X,1PE16.8,2X,I10,2X,1PE16.8,2X,1PE16.8)") LLbin, LLHisto(1)%LowVal+LLbin*LLHisto(1)%BinSize, LLHisto(1)%Hits(LLbin),LLHisto(2)%LowVal+LLbin*LLHisto(2)%BinSize, LLHisto(2)%Hits(LLbin),alpha(1,LLBin),alpha(2,LLBin)
   enddo


!************************************************************
!  5. Now find the point at which alpha=1-beta
!************************************************************

! for the time being, I'm going to assume that both distributions have the same range and binning - can change this later
   do LLbin=1,LLHisto(1)%NBins
      checkalpha(LLBin)=alpha(1,LLBin)+alpha(2,LLBin)-1d0
      !   print *, LLBIn,alpha(1,LLBin),alpha(2,LLBin),checkalpha(LLBin)
      if (checkalpha(LLBin)*checkalpha(LLBin-1) .lt. 0d0) then
         alphamin=alpha(1,LLBin-1)
         alphamax=alpha(1,LLBin)
         betamin =alpha(2,LLBin-1)
         betamax =alpha(2,LLBin)
         
         if (checkalpha(LLBin) .eq. 1d0 .and. checkalpha(LLBin-1) .eq.-1d0) then
            ! this is comparing two identical hypotheses!
            alphamin=0.5d0
            alphamax=0.5d0
            betamin =0.5d0
            betamax =0.5d0
         endif
      endif
   enddo

   
   alpha_ave=(alphamin+alphamax)/2d0
   beta_ave=(betamin+betamax)/2d0
   sigs=dsqrt(2d0)*inverf(1d0-alpha_ave)
   
   print *, 'alpha value in range:', alphamin,alphamax
   print *, 'beta value in range:', betamin,betamax 
   print *, 'alpha average :', alpha_ave
   print *, 'beta average :', beta_ave
   print *, 'sigma value :', sigs
   
   write(14,"(A,2X,1PE16.8,2X,1PE16.8)") "# alpha value in range:",alphamin,alphamax
   write(14,"(A,2X,1PE16.8,2X,1PE16.8)") "# beta value in range:",betamin,betamax
   write(14,"(A,2X,1PE16.8,2X,1PE16.8)") "# alpha average:", alpha_ave 
   write(14,"(A,2X,1PE16.8,2X,1PE16.8)") "# beta average:", beta_ave   
   write(14,"(A,2X,1PE16.8,2X,1PE16.8)") "# sigma value: ", sigs
   

   write(15,"(1PE16.8)") alpha_ave
   write(15,"(1PE16.8)") beta_ave
   write(15,"(1PE16.8)") sigs
   

endif! MPI_Rank

!DEC$ IF(_UseMPI .EQ.1)
   call MPI_FINALIZE(ierror)
!DEC$ ENDIF



contains

  FUNCTION Poisson(nu,n)
! Poisson distribution = exp(-nu)*nu^n/n!
    implicit none
    integer :: nu,n
    real(8) :: Poisson
    
    Poisson=-nu+n*log(1d0*nu)-logfac(n)
    Poisson=exp(Poisson)

  end FUNCTION POISSON


  FUNCTION logfac(N)
     ! log(N!)=log[ (N)(N-1)(N-2)...(2)(1)]=log(N)+log(N-1)+...+log(2)
    implicit none
    integer :: N,i
    real(8) :: logfac
    
    logfac=0d0
    do i=2,N
       logfac=logfac+dlog(1d0*i)
    enddo

  end FUNCTION logfac


  function INVERF(Y)
    real(8) :: Y,inverf,zork,inverf_old,inverf_new,derinverf_old,err
    integer :: Nmax,k
    real(8), parameter :: tol=1d-10
    real(8), parameter :: DblPi = 3.1415926535897932384626433832795028842d0
! done with Newton-Raphson method, using Maclaurin series as initial guess

! initial guess
    Nmax=10
    INVERF=0d0
    do k=0,Nmax
       zork=inverf_coeff(k)
       INVERF=INVERF + INVERF_COEFF(k)/(2*k+1)*( dsqrt(Dblpi)/2d0*Y )**(2*k+1)
    enddo
!    print *, 'first guess',inverf

! NR method
    inverf_old=inverf
    err=1d8
    k=0
    do while (err > tol .and. k .lt. 50)
       k=k+1
       inverf_old=inverf
       derinverf_old=2d0/dsqrt(DblPi)*exp(-inverf_old**2)
       inverf_new=inverf_old + (Y-erf(inverf_old))/derinverf_old
       err = abs(inverf_new-inverf_old)/inverf_old
       inverf=inverf_new
    enddo
    
    if (k .eq. 50) then
!       print *, 'INVERF REACHED MAXIMUM 50 NEWTON-RAPHSON ITERATIONS!'
    else
!       print *, 'NEWTON-RAPHSON CONVERGES WITH TOLERANCE ', tol, ' IN ', k, 'ITERATIONS'
    endif

!    print *, 'nr',inverf
    return
    

  end function INVERF       



  recursive function INVERF_COEFF(k)
    integer :: k,m
    real(8) :: inverf_coeff,c,inverf_koeff(0:44)
    
    c = 0d0;
    if (k .eq. 0) then
       inverf_coeff=1d0
       return
    endif
    do m=0,k-1
       c = c + INVERF_COEFF(m) *  INVERF_COEFF(k-1-m) / ( (m+1)*(2*m+1) )
    enddo
    INVERF_COEFF=c
    return

  end function INVERF_COEFF           


  subroutine init_random_seed()
! this ensures that each of the workers uses a random seed
! taken from http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid, t(2), s
    integer(8) :: count, tms
    
    call random_seed(size = n)
    allocate(seed(n))
    ! First try if the OS provides a random number generator
!    open(newunit=un, file="/dev/urandom", access="stream", &
!         form="unformatted", action="read", status="old", iostat=istat)
!    if (istat == 0) then
!       read(un) seed
!       close(un)
!    else
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       call system_clock(count)
       if (count /= 0) then
          t = transfer(count, t)
       else
          call date_and_time(values=dt)
          tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24 * 60 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
          t = transfer(tms, t)
       end if
       s = ieor(t(1), t(2))
       pid = getpid() + 1099279 ! Add a prime
       s = ieor(s, pid)
       if (n >= 3) then
          seed(1) = t(1) + 36269
          seed(2) = t(2) + 72551
          seed(3) = pid
          if (n > 3) then
             seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
          end if
       else
          seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
       end if
!    end if
    call random_seed(put=seed)
  end subroutine init_random_seed



end program BinnedLogL
  
