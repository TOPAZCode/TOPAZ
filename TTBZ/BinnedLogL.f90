program BinnedLogL
  use ifport
  implicit none
  integer :: NumArgs,Histo,NPseudoExp
  real(8) :: PreFactor,Data
  character :: operation*(10),H0_infile*(50),H1_infile*(50),Data_str*(5),Histo_str*(5),NPseudoExp_str*(9),dummy*(1),outfile*(50)
  character(len=*),parameter :: fmt1 = "(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)"
character(len=*),parameter :: fmt2 = "(I5,A,2X,1PE14.7,A,2X,I9,A,2X,1PE14.7,A,2X,I9,A,2X,1PE14.7,A,2X,1PE14.7,A,2X,1PE14.7,A)"
  integer,parameter :: MaxBins=1100, MaxEvents=100000
  integer,parameter :: NumLLbins=1100
  real(8) :: BinVal(1:2,1:MaxBins)=-1d-99,Value(1:2,1:MaxBins)=-1d-99,Error(1:2,1:MaxBins)=-1d-99,BinSize0,BinSize1
  integer :: NHisto(1:2)=-999999,Hits(1:2,1:MaxBins)=-999999,NumBins0=0,NumBins1=0,iPseudoExp,iBin,LLBin
  real(8) :: TotCS0,TotCS1,PoissonMax0(1:MaxBins),PoissonMax1(1:MaxBins)
  logical :: GotNumEvents(1:MaxBins)
  integer :: ExpectedEvents0(1:MaxBins),ExpectedEvents1(1:MaxBins)
  real(8) :: nran(1:2),offset
  integer :: TryEvts,ObsEvents0(1:MaxBins),ObsEvents1(1:MaxBins)
  integer :: PlotObsEvts(1:MaxBins,1:MaxBins)
  real(8) :: WhichBin,IntH0,IntH1,alpha(1:1000), beta(1:1000),check(1:1000),LLRatio
  type :: Histogram
     integer :: NBins
     real(8) :: BinSize
     real(8) :: LowVal
     integer :: Hits(1:1000)
  end type Histogram
  type(Histogram) :: LLHisto1,LLHisto2

  LLHisto1%NBins   = 1000
  LLHisto1%BinSize = 0.1d0
  LLHisto1%LowVal  = 0d0
  LLHisto1%Hits(:) = 0
  
  LLHisto2%NBins   = 1000
  LLHisto2%BinSize = 0.1d0
  LLHisto2%LowVal  = 0d0
  LLHisto2%Hits(:) = 0

  NumArgs=NArgs()-1
  PreFactor=8d0

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
  call GetArg(7,outfile)
  read(Histo_str,"(I2)") Histo
  read(Data_str,"(F10.6)") Data
  read(NPseudoExp_str,"(I9)") NPseudoExp

!--------------------------------------------------
!          1. reading input files
!--------------------------------------------------

  open(unit=12,file=trim(H0_infile),form='formatted',access='sequential')  ! open input file1
  open(unit=13,file=trim(H1_infile),form='formatted',access='sequential')  ! open input file1
  open(unit=14,file=trim(outfile),form='formatted',access='sequential')  ! open output file

  do while(.not.eof(12))!   reading input file 1
     read(unit=12,fmt="(A)") dummy
     if(dummy(1:1).eq."#") cycle
     backspace(unit=12) ! go to the beginning of the line
     
     read(unit=12,fmt=fmt1) NHisto(1),dummy,BinVal(1,NumBins0+1),dummy,Value(1,NumBins0+1),dummy,Error(1,NumBins0+1),dummy,Hits(1,NumBins0+1),dummy
     if( NHisto(1).ne.Histo ) cycle
     TotCS0=TotCS0+Value(1,NumBins0+1)*(BinVal(1,NumBins0+1)-BinVal(1,NumBins0))
     NumBins0=NumBins0 + 1
  enddo
  BinSize0=BinVal(1,NumBins0)-BinVal(1,NumBins0-1)
  
  do while(.not.eof(13))!   reading input file 2
     read(unit=13,fmt="(A)") dummy
     if(dummy(1:1).eq."#") cycle
     backspace(unit=13) ! go to the beginning of the line
     
     read(unit=13,fmt=fmt1) NHisto(2),dummy,BinVal(2,NumBins1+1),dummy,Value(2,NumBins1+1),dummy,Error(2,NumBins1+1),dummy,Hits(2,NumBins1+1),dummy
     if( NHisto(2).ne.Histo ) cycle
     TotCS1=TotCS1+Value(2,NumBins1+1)*(BinVal(2,NumBins1+1)-BinVal(2,NumBins1))
     NumBins1=NumBins1 + 1
     
  enddo
  BinSize1=BinVal(2,NumBins1)-BinVal(2,NumBins1-1)

  if( NumBins1.ne.NumBins0 ) then
     print *, "Error: Number of bins in input file 1 and 2 are different: ",NumBins0,NumBins1
     stop
  endif

! print the input histograms 
  write(*,"(2X,A,16X,A,11X,A,16X,A)") "NBin|","Input file 1","|","Input file 2"
  do iBin=1,NumBins0
    write(*,fmt="(2X,1I3,A,2X,1PE10.3,2X,1PE23.16,A,2X,1PE10.3,2X,1PE23.16)") iBin," | ",BinVal(1,iBin),Value(1,iBin)," | ",BinVal(2,iBin),Value(2,iBin)
    if( dabs(BinVal(1,iBin)-BinVal(2,iBin)).gt.1d-6 ) then
        print *, "Error: Different bin sizes in input files 1 and 2"
        stop
    endif
  enddo

! -----------------------------------------------------------------
! 2. Find the expected values for null and alt. hypothesis in each bin
! -----------------------------------------------------------------

! Null hypothesis first
  do iBin=1,NumBins0
     ExpectedEvents0(iBin)=int(Value(1,iBin)*Data*PreFactor*BinSize0)
     ExpectedEvents1(iBin)=int(Value(2,iBin)*Data*PreFactor*BinSize1)

     ! now find the value of the Poisson distribution at this maximum
     PoissonMax0(iBin)=Poisson(ExpectedEvents0(iBin),ExpectedEvents0(iBin))
     PoissonMax1(iBin)=Poisson(ExpectedEvents1(iBin),ExpectedEvents1(iBin))
  enddo
  
  call random_seed()

! -----------------------------------------------------------------
! 3.1 Generate Poisson distribution about expected null value in each bin
! -----------------------------------------------------------------

  PlotObsEvts=0
  do iPseudoExp=1,NPseudoExp
     GotNumEvents=.false.
     do iBin=1,NumBins0
        do while (.not. GotNumEvents(iBin))
           call random_number(nran(1:2))
           nran(1)=nran(1)*PoissonMax0(iBin)        
           TryEvts=int(3d0*ExpectedEvents0(iBin)*nran(2))
           
           if (Poisson(ExpectedEvents0(iBin),TryEvts) .gt. nran(1)) then
              ObsEvents0(iBin)=TryEvts
              GotNumEvents(iBin)=.true.
           endif
           if (ExpectedEvents0(iBin) .eq. 0) then
              ObsEvents0(iBin)=0
              GotNumEvents(iBin)=.true.
           endif
        enddo
        PlotObsEvts(iBin,ObsEvents0(iBin))=PlotObsEvts(iBin,ObsEvents0(iBin))+1
     enddo

! -----------------------------------------------------------------
! 3.2 Now find the log likelihood, with a Poisson distr in each bin
! -----------------------------------------------------------------

     LLRatio=0d0
     do iBin=1,NumBins0
        if (ObsEvents0(iBin) .ne. 0 .and. ExpectedEvents0(iBin) .ne. 0 .and. &
             & ExpectedEvents1(iBin) .ne. 0) then
           LLRatio=LLRatio+ObsEvents0(iBin)*dlog(1d0*ExpectedEvents0(iBin)/ExpectedEvents1(iBin))
           
        endif
     enddo

! offset to get the distributions with the histogram limits
     if (iPseudoExp .eq. 1) then
        offset=-100d0*(int(LLRatio)/100)
     endif
     LLRatio=LLRatio+offset

     if (iPseudoExp .eq. 1) then
        if (LLRatio .lt. LLHisto1%LowVal .or. LLRatio .gt. (LLHisto1%LowVal+LLHisto1%NBins*LLHisto1%BinSize)) then
            print *, 'LLRatio out  of binning range! Reset binning and rerun'
            print *, LLRatio
            stop
         endif
      endif

! -----------------------------------------------------------------
! 3.3 Bin the log likelihood
! -----------------------------------------------------------------

      WhichBin = (LLRatio-LLHisto1%LowVal)/LLHisto1%BinSize + 1
      if( WhichBin.lt.0 ) WhichBin = 1
      if( WhichBin.gt.LLHisto1%NBins ) WhichBin = LLHisto1%NBins
      LLHisto1%Hits(WhichBin) = LLHisto1%Hits(WhichBin) + 1
      

   enddo! iPseudoExp


   do iPseudoExp=1,1000    ! this is just a check 
      write(201,*) iPseudoExp,PlotObsEvts(1,iPseudoExp)
   enddo

! -----------------------------------------------------------------
! 4. Calculate the integral of the log likelihood ratio curve, at each bin
! -----------------------------------------------------------------

  do LLbin=1,LLHisto1%NBins
     IntH0=IntH0+LLHisto1%BinSize*LLHisto1%Hits(LLbin)
     alpha(LLbin)=IntH0/(NPseudoExp*LLHisto1%BinSize)
     print *, LLbin, LLHisto1%LowVal + LLbin*LLHisto1%BinSize, LLHisto1%Hits(LLbin), alpha(LLBin)
  enddo



! Now alt. hypothesis

! -----------------------------------------------------------------
! 5.1 Generate Poisson distribution about expected null value in each bin
! -----------------------------------------------------------------

  PlotObsEvts=0
  do iPseudoExp=1,NPseudoExp
     GotNumEvents=.false.
     do iBin=1,NumBins1
        do while (.not. GotNumEvents(iBin))
           call random_number(nran(1:2))
           nran(1)=nran(1)*PoissonMax1(iBin)        
           TryEvts=int(3d0*ExpectedEvents1(iBin)*nran(2))

           if (Poisson(ExpectedEvents1(iBin),TryEvts) .gt. nran(1)) then
              ObsEvents1(iBin)=TryEvts
              GotNumEvents(iBin)=.true.
           endif
           if (ExpectedEvents1(iBin) .eq. 0) then
              ObsEvents1(iBin)=0
              GotNumEvents(iBin)=.true.
           endif
        enddo
        PlotObsEvts(iBin,ObsEvents1(iBin))=PlotObsEvts(iBin,ObsEvents1(iBin))+1
     enddo

! -----------------------------------------------------------------
! 5.2 Now find the log likelihood, with a Poisson distr in each bin
! -----------------------------------------------------------------
   
     LLRatio=0d0
     do iBin=1,NumBins1
        if (ObsEvents1(iBin) .ne. 0 .and. ExpectedEvents0(iBin) .ne. 0 .and. &
             & ExpectedEvents1(iBin) .ne. 0) then
           LLRatio=LLRatio+ObsEvents1(iBin)*dlog(1d0*ExpectedEvents0(iBin)/ExpectedEvents1(iBin))
        endif
     enddo
     LLRatio=LLRatio+offset

     if (iPseudoExp .eq. 1) then
        if (LLRatio .lt. LLHisto2%LowVal .or. LLRatio .gt. (LLHisto2%LowVal+LLHisto2%NBins*LLHisto2%BinSize)) then
            print *, 'LLRatio out  of binning range! Reset binning and rerun'
            print *, LLRatio
            stop
         endif
      endif

! -----------------------------------------------------------------
! 5.3 Bin the log likelihood
! -----------------------------------------------------------------

      WhichBin = (LLRatio-LLHisto2%LowVal)/LLHisto2%BinSize + 1
      if( WhichBin.lt.0 ) WhichBin = 1
      if( WhichBin.gt.LLHisto2%NBins ) WhichBin = LLHisto2%NBins
      LLHisto2%Hits(WhichBin) = LLHisto2%Hits(WhichBin) + 1

   enddo! iPseudoExp


   do iPseudoExp=1,1000     ! again, a check
      write(202,*) iPseudoExp,PlotObsEvts(1,iPseudoExp)
   enddo

! -----------------------------------------------------------------
! 6. Calculate the integral of the log likelihood ratio curve, at each bin, and !    write out
! -----------------------------------------------------------------

   do LLbin=1,LLHisto2%NBins
      IntH1=IntH1+LLHisto2%BinSize*LLHisto2%Hits(LLbin)
      beta(LLbin)=IntH1/(NPseudoExp*LLHisto2%BinSize)
      print *, LLbin, LLHisto2%LowVal + LLbin*LLHisto2%BinSize, LLHisto2%Hits(LLbin),beta(LLBin)
      check(LLbin)=alpha(LLbin)+beta(LLbin)-1d0
      
      if (LLHisto1%Hits(LLbin) .ne. 0 .or. LLHisto2%Hits(LLbin) .ne. 0) then
         write(14,fmt2) LLbin,"|", LLHisto1%LowVal + LLbin*LLHisto1%BinSize, "|", LLHisto1%Hits(LLbin), "|", LLHisto2%LowVal + LLbin*LLHisto2%BinSize, "|", LLHisto2%Hits(LLbin), "|", alpha(LLBin), "|", beta(LLbin), "|",alpha(LLbin)+beta(LLbin)-1d0, "|"
      endif
   enddo
  
! -----------------------------------------------------------------
! 7. Find the point at which alpha=1-beta -- this is the p value
! ----------------------------------------------------------------- 

   do LLBin=1,LLHisto2%NBins
      if (check(LLBin)*check(LLBin+1) .le. 0d0) then
         print *, alpha(LLBin), beta(LLBin)
         print *, alpha(LLBin+1), beta(LLBin+1)
         stop
      endif
   enddo

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
       
           
end program BinnedLogL
  
