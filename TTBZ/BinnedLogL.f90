! ./BinnedLogL GetLogL LO_Hi17_V+0.00_A+0.00.dat LO_Hi17_V+4.00_A+0.20.new.dat 17 3000 20000 1.2 test
! Compile:  ifort -o BinnedLogL -D_UseMPI=0 BinnedLogL.f90
! Compile:  mpif90 -D_UseMPI=1 -o BinnedLogL BinnedLogL.f90
! Compile:  safety options:   
  PROGRAM BinnedLogL
  use ifport
  implicit none
  integer :: NumArgs,Histo,NPseudoExp,NPseudoExp_Worker,Data
  real(8) :: PreFactor
  character :: operation*(10),H0_infile*(100),H1_infile*(100),Data_str*(5),Histo_str*(5),NPseudoExp_str*(9),dummy*(1),outfile*(100),DeltaN_str*(5),filename_out1*(100),filename_out2*(100)
  character(len=*),parameter :: fmt1 = "(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)"
  character(len=*),parameter :: fmt2 = "(I5,A,2X,1PE14.7,A,2X,I9,A,2X,1PE14.7,A,2X,I9,A,2X,1PE14.7,A,2X,1PE14.7,A,2X,1PE14.7,A)"
  integer,parameter :: MaxBins=100, MaxEvents=100000,MaxExp=1000000
  real(8) :: BinVal(1:2,1:MaxBins)=-1d-99,Value(1:2,1:MaxBins)=-1d-99,Error(1:2,1:MaxBins)=-1d-99,BinSize(1:2)
  integer :: NHisto(1:2)=-999999,Hits(1:2,1:MaxBins)=-999999,NumBins0=0,NumBins1=0,NumBins,BinMin,BinMax,iPseudoExp,iBin,iHypothesis,LLBin
  logical :: GotNumEvents(1:MaxBins),IdenticalHypothesis
  integer :: ExpectedEvents(1:2,1:MaxBins),TryEvts,MinEvts,MaxEvts,ObsEvents(1:2,1:MaxBins),PlotObsEvts(1:2,1:10000)=0d0,TotSigmaHisto
  real(8) :: nran(1:2),offset,sran,DeltaN,Delta(1:2),MaxAllowedDelta,AvgTotalEvents,LLRatio,PoissonWidth
  real(8) :: alphamin,alphamax,betamin,betamax,rescale(1:2),alpha_ave,beta_ave,sigs
  real(8) :: LLRatio_array(1:2,1:MaxExp),LLRatio_array_tmp(1:2,1:MaxExp),LLRatio_min,LLRatio_max,WhichBin
  integer :: j,i,s,MPI_Rank,worker,MPI_NUM_PROCS,ierror,TotalEvents(1:2)=0,TheUnit
  integer,parameter :: MaxLLBins=1000
  real(8) :: sigmatot(1:2),check(1:2),alpha(1:MaxLLBins),beta(1:MaxLLBins),IntLLRatio(1:2),checkalpha(0:MaxLLBins)=0d0
  type :: Histogram
     integer :: NBins
     real(8) :: BinSize
     real(8) :: LowVal
     integer :: Hits(1:MaxLLBins)
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



  NumArgs=NArgs()-1
  PreFactor=8d0!  Lepton multiplicities


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
  read(Data_str,"(I5)") Data
  read(NPseudoExp_str,"(I9)") NPseudoExp
  read(DeltaN_str,"(F5.1)") DeltaN
  if (NPseudoExp .gt. MaxExp) then
     print *, "Too many experiments!"
     print *, MaxExp, NPseudoExp
     stop
  endif

  if( DeltaN.lt.0d0 ) then
     print *, "!! WARNING DeltaN has to be positive",DeltaN
     stop
  endif
  if( DeltaN.gt.100.0d0) then
     print *, "!! WARNING DeltaN is outside the range 0..100%",DeltaN
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

  if( MPI_Rank.eq.0 ) then
  do TheUnit=6,14,8! write to screen(unit=6) and test.dat(unit=14)
    write(TheUnit,"(A)") "# Performing binned log-likelihood analysis."
    write(TheUnit,"(A,I2,A,I5,A,I9,A,F6.1,A)") "# Reading histogram ",Histo," for likelihood analysis with Lumi=",Data," fb^-1 using ",NPseudoExp," pseudo-experiments, and scale uncertainty", DeltaN,"%."
    write(TheUnit,"(A,A)") '# Input file 1: NULL hypothesis        ', trim(H0_infile)
    write(TheUnit,"(A,A)") '# Input file 2: ALTERNATIVE hypothesis ', trim(H1_infile)
    write(TheUnit,"(A,A,A,A)") "# Writing output to files : ", trim(filename_out1), " and ", trim(filename_out2)
    write(TheUnit,"(A)") "#";  write(TheUnit,"(A)") "#";
  enddo
  endif

!-------------------------------------------------------------------------
! 1a. Because our Delta phi_ll histogram excludes a small bin at the end,
! first use histo #TotSigmaHisto  to get the total cross-section
!-------------------------------------------------------------------------

  TotSigmaHisto = 1 ! Histo

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
    if( NHisto(2).ne.TotSigmaHisto ) cycle
     NumBins1=NumBins1 + 1
  enddo

  BinSize(1) = BinVal(1,2)-BinVal(1,1)
  BinSize(2) = BinVal(2,2)-BinVal(2,1)
  do iBin=1,NumBins1! calculate total cross section
       sigmatot(1) = sigmatot(1) + Value(1,iBin) * BinSize(1)
       sigmatot(2) = sigmatot(2) + Value(2,iBin) * BinSize(2)
  enddo

  if( MPI_Rank.eq.0 ) then
  do TheUnit=6,14,8! write to screen(unit=6) and test.dat(unit=14)
      write(TheUnit,"(A,1PE16.8,A,I6)") "# Total cross section of input file 1: ",sigmatot(1),"   <-->   Number of events: ",int(Data * sigmatot(1) * PreFactor)
      write(TheUnit,"(A,1PE16.8,A,I6)") "# Total cross section of input file 2: ",sigmatot(2),"   <-->   Number of events: ",int(Data * sigmatot(2) * PreFactor) 
  enddo
  endif



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
 

! select the range of used histogram binsize  
  BinMin=1            +1 ! remove first and last bin for NLO analysis
  BinMax=NumBins      -1
  if( BinMin.ne.1 .or. BinMax.ne.NumBins ) then
      do TheUnit=6,14,8! write to screen(unit=6) and test.dat(unit=14)
        if( MPI_Rank.eq.0 ) write(TheUnit,"(A,I2,A,I2)") "# WARNING: Not all histogram bins are used for the analysis. BinMin,BinMax=",BinMin,",",BinMax
      enddo
  endif

! print the input histograms 
  if( MPI_Rank.eq.0 ) write(*,"(A,16X,A,11X,A,16X,A)") "# NBin|","Input file 1","|","Input file 2"
  do iBin=BinMin,BinMax
     if( MPI_Rank.eq.0 ) write(*,fmt="(2X,1I3,A,2X,1PE10.3,2X,1PE23.16,A,2X,1PE10.3,2X,1PE23.16)") iBin," | ",BinVal(1,iBin),Value(1,iBin)," | ",BinVal(2,iBin),Value(2,iBin)
     if( dabs(BinVal(1,iBin)-BinVal(2,iBin)).gt.1d-6 ) then
        print *, "Error: Different bin sizes in input files 1 and 2"
        stop
     endif
  enddo

  
! -----------------------------------------------------------------
! 2. Find the expected values for null and alt. hypothesis in each bin
! -----------------------------------------------------------------

  do iHypothesis=1,2    ! 1=null, 2=alt
     TotalEvents(iHypothesis) = int( sum(Value(iHypothesis,BinMin:BinMax))*Data*PreFactor*BinSize(iHypothesis) )
     if( MPI_Rank.eq.0 ) write(*,"(A,I1,A,I6)") "Hypothesis:",iHypothesis," Number of events in Histogram:",TotalEvents(iHypothesis)
  enddo
  if( ( TotalEvents(1).lt.100) .or. ( TotalEvents(2).lt.100) ) print *, "!! WARNING: less than 100 events in total histogram !!"
 
  if( TotalEvents(1).gt.TotalEvents(2) ) then! null hypothesis has larger total cross section
        MaxAllowedDelta = 0.5d0*(dble(TotalEvents(1))/dble(TotalEvents(2))-1d0) * 100d0
        AvgTotalEvents  = 0.5d0*(TotalEvents(1)+TotalEvents(2))
        if( DeltaN .lt. MaxAllowedDelta ) then! if DeltaN is smaller than difference between the two values
            Delta(1) = -DeltaN/100d0! make null hypothesis smaller
            Delta(2) = +DeltaN/100d0! make alt. hypothesis larger
        else
            Delta(1) = dble(AvgTotalEvents)/dble(TotalEvents(1)) - 1d0! rescale such that total cross section is equal to the average between the two hypothesis
            Delta(2) = dble(AvgTotalEvents)/dble(TotalEvents(2)) - 1d0
        endif
  else! null hypothesis has smaller total cross section
        MaxAllowedDelta = 0.5d0*(dble(TotalEvents(2))/dble(TotalEvents(1))-1d0) * 100d0
        AvgTotalEvents  = 0.5d0*(TotalEvents(1)+TotalEvents(2))
        Delta(1) = +DeltaN/100d0! make null hypothesis larger
        Delta(2) = -DeltaN/100d0! make alt. hypothesis smaller
        if( dble(TotalEvents(1))*(1d0+Delta(1)) .gt.  dble(TotalEvents(2))*(1d0+Delta(2)) ) then
            Delta(1) = dble(AvgTotalEvents)/dble(TotalEvents(1)) - 1d0! rescale such that total cross section is equal to the average between the two hypothesis
            Delta(2) = dble(AvgTotalEvents)/dble(TotalEvents(2)) - 1d0
        endif
  endif
  Value(1,:) = Value(1,:) * ( 1d0 + Delta(1) )
  Value(2,:) = Value(2,:) * ( 1d0 + Delta(2) )
  do iHypothesis=1,2    ! 1=null, 2=alt
     TotalEvents(iHypothesis) = int( sum(Value(iHypothesis,BinMin:BinMax))*Data*PreFactor*BinSize(iHypothesis) )
     if( MPI_Rank.eq.0 ) write(*,"(A,I1,A,I6,A,F5.1,A)") "Hypothesis:",iHypothesis," Number of events in Histogram:",TotalEvents(iHypothesis)," (after rescaling by uncertainty of ",Delta(iHypothesis)*100d0,"%)"
  enddo


  
!  call random_seed()
  call init_random_seed()
! -----------------------------------------------------------------
! 3.1 Generate Poisson distribution about expected null value in each bin
! -----------------------------------------------------------------
  NPseudoExp_Worker = NPseudoExp/(MPI_NUM_PROCS)

!DEC$ IF(_UseMPI .EQ.1)
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
!DEC$ ENDIF
!  call random_seed()
  do iHypothesis=1,2
     write(*,"(A,I2,A,I2,A,I7,A)") 'Hypothesis=', iHypothesis, ' on worker=',MPI_Rank," running ",NPseudoExp_Worker," pseudoexperiments" 
     do iPseudoExp=1,NPseudoExp_Worker
!        if( mod(iPseudoExp,25000).eq.0 ) print *, "Pseudo experiment ",iPseudoExp,"/",NPseudoExp

        ExpectedEvents(1,:)=int(Value(1,:)*Data*PreFactor*BinSize(1))
        ExpectedEvents(2,:)=int(Value(2,:)*Data*PreFactor*BinSize(2))

        GotNumEvents=.false.
        do iBin=BinMin,BinMax
           do while (.not. GotNumEvents(iBin))
              call random_number(nran(1:2))
              nran(1) = nran(1) * Poisson(ExpectedEvents(iHypothesis,iBin),ExpectedEvents(iHypothesis,iBin))
!               MinEvts=0
!               MaxEvts=5*ExpectedEvents(iHypothesis,iBin)! old: linear distribution between 0 and 5*(Poisson maximum)
!               TryEvts=int((MaxEvts-MinEvts)*nran(2)+MinEvts)             
              MinEvts = max(0d0,ExpectedEvents(iHypothesis,iBin) - 3.5d0*dsqrt(dble(ExpectedEvents(iHypothesis,iBin)))) ! generate events between +/-3.5 sigma around Poisson maximum
              MaxEvts = ExpectedEvents(iHypothesis,iBin)         + 3.5d0*dsqrt(dble(ExpectedEvents(iHypothesis,iBin)))
              TryEvts=int((MaxEvts-MinEvts)*nran(2)+MinEvts)
              if( Poisson(ExpectedEvents(iHypothesis,iBin),TryEvts) .gt. nran(1) ) then
                 ObsEvents(iHypothesis,iBin)=TryEvts
                 GotNumEvents(iBin)=.true.
              endif
              if (ExpectedEvents(iHypothesis,iBin) .eq. 0) then
                 ObsEvents(iHypothesis,iBin)=0
                 GotNumEvents(iBin)=.true.
              endif
           enddo

!            if( iBin.eq.1 .and. ObsEvents(iHypothesis,iBin).gt.0 ) then! this is to check the Poisson distribution in each bin
!               PlotObsEvts(iHypothesis,ObsEvents(iHypothesis,iBin))=PlotObsEvts(iHypothesis,ObsEvents(iHypothesis,iBin))+1
!            endif
        enddo


! -----------------------------------------------------------------
! 3.2 Now find the log likelihood, with a Poisson distr in each bin
! -----------------------------------------------------------------
     LLRatio=0d0
     do iBin=BinMin,BinMax
        if (ObsEvents(iHypothesis,iBin) .ne. 0 .and. ExpectedEvents(1,iBin) .ne. 0 .and. ExpectedEvents(2,iBin) .ne. 0) then
              LLRatio = LLRatio + ObsEvents(iHypothesis,iBin)*dlog(dble(ExpectedEvents(1,iBin))/dble(ExpectedEvents(2,iBin))) - dble(ExpectedEvents(1,iBin)) + dble(ExpectedEvents(2,iBin))
        endif
     enddo

     
! offset to get the distributions with the histogram limits
!DEC$ IF(_UseMPI .NE.1)
     if (iHypothesis .eq. 1 .and. iPseudoExp .eq. 1) then
        write(*,*) 'LLRatio(before offset)=',LLRatio
        offset=-100d0*(int(LLRatio)/100)
        write(*,*) 'using offset of ', offset
     endif  
     LLRatio=LLRatio+offset
     if (iPseudoExp .eq. 1) then
        write(*,*) 'LLRatio(after offset)=',LLRatio
     endif
!DEC$ ENDIF
     LLRatio_array(iHypothesis,iPseudoExp)=LLRatio
        
  enddo!  iPseudoExp

!   do i=1,1000
!      s=101+iHypothesis
!      write(s,*) i,PlotObsEvts(iHypothesis,i),ExpectedEvents(iHypothesis,1)
!   enddo

enddo!   iHypothesis


  write(*,*) "MPI rank ",MPI_Rank," finished generating pseudoexperiments"



!DEC$ IF(_UseMPI .EQ.1)
   call MPI_BARRIER(MPI_COMM_WORLD,ierror)
   if( MPI_Rank.eq.0 ) then
      do worker=1,MPI_NUM_PROCS-1
        call MPI_Recv(LLRatio_array_tmp, 2*MaxExp,MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE ,1,MPI_COMM_WORLD, status,ierror)
        write(*,*) "received array from worker",worker

        LLRatio_array(1,NPseudoExp_Worker*worker+1:NPseudoExp_Worker*(worker+1)) = LLRatio_array_tmp(1,1:NPseudoExp_Worker)
        LLRatio_array(2,NPseudoExp_Worker*worker+1:NPseudoExp_Worker*(worker+1)) = LLRatio_array_tmp(2,1:NPseudoExp_Worker)
        if (worker.eq.1) then
           offset=-100d0*(int(LLRatio_array(1,1))/100)
           write(*,*) "using offset of ", offset
        endif
      enddo
      LLRatio_array=LLRatio_array+offset

   elseif( MPI_Rank.gt.0 ) then
        call MPI_Send(LLRatio_array, 2*MaxExp,MPI_DOUBLE_PRECISION, 0 ,1,MPI_COMM_WORLD,ierror)
   endif
!DEC$ ENDIF







if( MPI_Rank.eq.0 ) then 

! -----------------------------------------------------------------
! 4.1 find a good bin size and range for histograms
! -----------------------------------------------------------------
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
   
   LLHisto(1)%NBins   = MaxLLBins
   LLHisto(1)%BinSize = (LLRatio_max-LLRatio_min)/LLHisto(1)%NBins
   LLHisto(1)%LowVal  = LLRatio_min
   LLHisto(1)%Hits(:) = 0
   LLHisto(2)%NBins   = LLHisto(1)%NBins
   LLHisto(2)%BinSize = LLHisto(1)%BinSize
   LLHisto(2)%LowVal  = LLHisto(1)%LowVal
   LLHisto(2)%Hits(:) = LLHisto(1)%Hits(:)
!    Write(*,*) 'using binsize', LLHisto(1)%BinSize, 'within',LLRatio_min, LLRatio_max

  
! -----------------------------------------------------------------
! 4.1 do the binning
! -----------------------------------------------------------------
   IdenticalHypothesis=.false.
   do iHypothesis=1,2
      do iPseudoExp=1,NPseudoExp
         LLRatio=LLRatio_array(iHypothesis,iPseudoExp)
         WhichBin = (LLRatio-LLHisto(iHypothesis)%LowVal)/LLHisto(iHypothesis)%BinSize + 1
         WhichBin=int(WhichBin)
         if( WhichBin.lt.0 ) WhichBin = 1
         if( WhichBin.gt.LLHisto(iHypothesis)%NBins ) WhichBin = LLHisto(iHypothesis)%NBins
         LLHisto(iHypothesis)%Hits(WhichBin) = LLHisto(iHypothesis)%Hits(WhichBin) + 1
      enddo

      if( any( LLHisto(iHypothesis)%Hits(:).eq.NPseudoExp ) ) then
           IdenticalHypothesis=.true.
      endif
   enddo



! -----------------------------------------------------------------
! 5. Calculate the integral of the log likelihood ratio curve, at each bin
! -----------------------------------------------------------------
!       do LLbin=1,LLHisto(iHypothesis)%NBins
!          IntLLRatio(iHypothesis)=IntLLRatio(iHypothesis)  +  LLHisto(iHypothesis)%BinSize * LLHisto(iHypothesis)%Hits(LLbin)
!          alpha(iHypothesis,LLBin)=IntLLRatio(iHypothesis)/(NPseudoExp*LLHisto(iHypothesis)%BinSize)
!       enddo
!    enddo
!    MARKUS: check that beta=1-alpha, check that intllratio/(NPseudoExp*LLHisto(iHypothesis)%BinSize) is one.

!  do LLbin=1,LLHisto1%NBins
!     IntH0=IntH0+LLHisto1%BinSize*LLHisto1%Hits(LLbin)
!     alpha(LLbin)=IntH0/(NPseudoExp*LLHisto1%BinSize)
!     print *, LLbin, LLHisto1%LowVal + LLbin*LLHisto1%BinSize, LLHisto1%Hits(LLbin), alpha(LLBin)
!  enddo




      iHypothesis=1
      IntLLRatio(:)=0
      do LLbin=1,LLHisto(iHypothesis)%NBins
         IntLLRatio(iHypothesis)=IntLLRatio(iHypothesis)  + LLHisto(iHypothesis)%BinSize * LLHisto(iHypothesis)%Hits(LLbin)
         alpha(LLBin)=IntLLRatio(iHypothesis)/(NPseudoExp*LLHisto(iHypothesis)%BinSize)! normalize total integral to one
      enddo

      iHypothesis=2
      IntLLRatio(:)=0
      do LLbin=LLHisto(iHypothesis)%NBins,1,-1
         IntLLRatio(iHypothesis)=IntLLRatio(iHypothesis)  + LLHisto(iHypothesis)%BinSize * LLHisto(iHypothesis)%Hits(LLbin)
         beta(LLBin)=IntLLRatio(iHypothesis)/(NPseudoExp*LLHisto(iHypothesis)%BinSize)! normalize total integral to one
      enddo


      write(*,*) "Writing log-likelihood distribution to output file"
      do LLbin=1,LLHisto(1)%NBins
          write(14,"(2X,I4,2X,1PE16.8,2X,I10,2X,1PE16.8,2X,I10,2X,1PE16.8,2X,1PE16.8)") LLbin, LLHisto(1)%LowVal+LLbin*LLHisto(1)%BinSize, LLHisto(1)%Hits(LLbin), &
                                                                                               LLHisto(2)%LowVal+LLbin*LLHisto(2)%BinSize, LLHisto(2)%Hits(LLbin), &
                                                                                               alpha(LLBin),beta(LLBin)
      enddo



!************************************************************
!  6. Now find the point at which alpha=beta
!************************************************************

! ! for the time being, I'm going to assume that both distributions have the same range and binning - can change this later
!    do LLbin=1,LLHisto(1)%NBins
!       checkalpha(LLBin)=alpha(1,LLBin)+alpha(2,LLBin)-1d0
!       !   print *, LLBIn,alpha(1,LLBin),alpha(2,LLBin),checkalpha(LLBin)
!       if (checkalpha(LLBin)*checkalpha(LLBin-1) .lt. 0d0) then
!          alphamin=alpha(1,LLBin-1)
!          alphamax=alpha(1,LLBin)
!          betamin =alpha(2,LLBin-1)
!          betamax =alpha(2,LLBin)
!          
!          if (checkalpha(LLBin) .eq. 1d0 .and. checkalpha(LLBin-1) .eq.-1d0) then
!             ! this is comparing two identical hypotheses!
!             alphamin=0.5d0
!             alphamax=0.5d0
!             betamin =0.5d0
!             betamax =0.5d0
!          endif
!       endif
!    enddo


      if( .not. IdenticalHypothesis ) then
        do LLbin=1,LLHisto(1)%NBins
          checkalpha(LLBin)=alpha(LLBin)-beta(LLBin)
          if( checkalpha(LLBin)*checkalpha(LLBin-1) .lt. 0d0 ) then! search for sign change
            alphamin=alpha(LLBin-1)
            alphamax=alpha(LLBin)
            betamin =beta(LLBin-1)
            betamax =beta(LLBin)
            exit! stop searching for sign change
          endif
        enddo
      else
        alphamin=0.5d0
        alphamax=0.5d0
        betamin =0.5d0
        betamax =0.5d0           
      endif

   
      alpha_ave=(alphamin+alphamax)/2d0
      beta_ave=(betamin+betamax)/2d0
      sigs=dsqrt(2d0)*inverf(1d0-alpha_ave)!  check this
      
      do TheUnit=6,14,8! write to screen(unit=6) and test.dat(unit=14)
        if( IdenticalHypothesis ) write(TheUnit,"(A)") "# NOTE: Identical Hypothesis!"
        write(TheUnit,"(A,2X,1PE16.8,2X,1PE16.8)") "# alpha value in range:",alphamin,alphamax
        write(TheUnit,"(A,2X,1PE16.8,2X,1PE16.8)") "# beta value in range:",betamin,betamax
        write(TheUnit,"(A,2X,1PE16.8,2X,1PE16.8)") "# alpha average:", alpha_ave 
        write(TheUnit,"(A,2X,1PE16.8,2X,1PE16.8)") "# beta average:", beta_ave   
        write(TheUnit,"(A,2X,1PE16.8,2X,1PE16.8)") "# sigma value: ", sigs
      enddo
      write(15,"(1PE16.8)") alpha_ave
      write(15,"(1PE16.8)") beta_ave
      write(15,"(1PE16.8)") sigs
   
      close(14)
      close(15)
      write(*,*) "Done"; write(*,*) "";  write(*,*) ""
endif! MPI_Rank



!DEC$ IF(_UseMPI .EQ.1)
   call MPI_FINALIZE(ierror)
!DEC$ ENDIF






!---------------------------------------------------------------------------------------------
CONTAINS
!---------------------------------------------------------------------------------------------











  FUNCTION Poisson(nu,n)
! Poisson distribution = exp(-nu)*nu^n/n!
    implicit none
    integer :: nu,n
    real(8) :: Poisson
    
    
    if( nu.le.0 ) then
        Poisson=0d0
    else
        Poisson=-nu+n*log(1d0*nu)-logfac(n)
        Poisson=exp(Poisson)
    endif

  RETURN
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



END PROGRAM
  
