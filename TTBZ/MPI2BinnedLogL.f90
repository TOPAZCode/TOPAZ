program BinnedLogL
  use ifport
  use mpi
  implicit none
  integer :: NumArgs,Histo,NPseudoExp
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
  real(8) :: LLRatio_array(1:2,1:MaxExp),LLRatio_min,LLRatio_max,WhichBin
    integer :: j,i,s,SUA
  type :: Histogram
     integer :: NBins
     real(8) :: BinSize
     real(8) :: LowVal
     integer :: Hits(1:1000)
  end type Histogram
  type(Histogram) :: LLHisto(1:2)
  ! this is all mpi stuff:
  integer :: my_id, root_process, ierr, status(MPI_STATUS_SIZE)
  integer :: num_procs,slaves,islave
  integer :: MA_InitPseudoExp,MA_FinalPseudoExp,MA_NPseudoExp,inittag,finaltag
  integer :: MAfSL_InitPseudoExp,MAfSL_FinalPseudoExp,MAfSL_NPseudoExp
  integer :: SL_InitPseudoExp,SL_FinalPseudoExp,SL_NPseudoExp,LLtag1,LLtag2,SL_inittag,SL_finaltag
  integer :: ThisSlaveInitPseudoExp,ThisSlaveFinalPseudoExp,ThisSlaveNPseudoExp,SlavesNPseudoExp
  real(8) :: LL1fromSlaves(1:MaxExp),LL2fromSlaves(1:MaxExp), SL_LLRatio_array(1:2,1:MaxExp),SL_LL2Ratio_array(1:MaxExp),SL_LL1Ratio_array(1:MaxExp)
  real(8) :: LL1_int(1:MaxExp),MA_LL1_int(1:MaxExp)



  call MPI_INIT (ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)
  root_process=0
  slaves=num_procs-1


  if (my_id .eq. root_process) then   ! do all the reading in, etc
  
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
        
        if( NHisto(1).ne.1 ) cycle
        NumBins0=NumBins0 + 1
     enddo
     
     do while(.not.eof(113))!   reading input file 2
        read(unit=113,fmt="(A)") dummy
        if(dummy(1:1).eq."#") cycle
        backspace(unit=113) ! go to the beginning of the line
        
        read(unit=113,fmt=fmt1) NHisto(2),dummy,BinVal(2,NumBins1+1),dummy,Value(2,NumBins1+1),dummy,Error(2,NumBins1+1),dummy,Hits(2,NumBins1+1),dummy
        if( NHisto(2).ne.1 ) cycle
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
     ! 1b. Now read in the actual ditribution of interest
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
           ExpectedEvents(iHypothesis,iBin)=int(Value(iHypothesis,iBin)*Data*PreFactor*BinSize(iHypothesis))
           if (ExpectedEvents(iHypothesis,iBin).eq.0) then
              PoissonMax(iHypothesis,iBin)=0d0
           else
              PoissonMax(iHypothesis,iBin)=Poisson(ExpectedEvents(iHypothesis,iBin),ExpectedEvents(iHypothesis,iBin))
           endif
        enddo
     enddo

     ! -----------------------------------------------------------------
     ! MPI: Assign each slave a small number of pseudoexpts
     ! -----------------------------------------------------------------
     
        print *, 'USING ', slaves, ' SLAVES' 
     do islave=1,slaves
        MA_NPseudoExp=NPseudoExp/slaves
        MA_InitPseudoExp=1+MA_NPseudoExp*(islave-1)
        MA_FinalPseudoExp=MA_InitPseudoExp+MA_NPseudoExp-1
        if (islave .eq. slaves) then   ! the last slave has to do more work!
           MA_FinalPseudoExp=NPseudoExp
           MA_NPseudoExp=MA_FinalPseudoExp-MA_InitPseudoExp+1
        endif
        write(110,*) , 'slave #: ', islave, 'jobs: ', MA_InitPseudoExp,MA_FinalPseudoExp,MA_NPseudoExp
        inittag=999
        finaltag=998
        call MPI_SEND(MA_InitPseudoExp,1,MPI_INT,islave,inittag,MPI_COMM_WORLD,ierr)
        call MPI_SEND(MA_NPseudoExp,1,MPI_INT,islave,inittag,MPI_COMM_WORLD,ierr)
     enddo
     

     do islave=1,slaves     
!        call MPI_RECV(MAfSL_InitPseudoExp,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
!        call MPI_RECV(MAfSL_InitPseudoExp,1,MPI_INT,islave,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
        write(110,*), 'receiving from slave #', islave
!        call MPI_RECV(MAfSL_NPseudoExp,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
!        call MPI_RECV(MAfSL_NPseudoExp,1,MPI_INT,islave,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
        MAfSL_FinalPseudoExp=MAfSL_InitPseudoExp+MAfSL_NPseudoExp-1
        LL1fromSlaves=1d0
!        call MPI_RECV(LL1fromSlaves(1),MAfSL_NPseudoExp,MPI_REAL,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
!        call MPI_RECV(LL2fromSlaves(1),MAfSL_NPseudoExp,MPI_REAL,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)

        call MPI_RECV(LL1fromSlaves(1),MAfSL_NPseudoExp,MPI_REAL,islave,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
!        call MPI_RECV(LL1fromSlaves(1),8,MPI_REAL,islave,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
!       call MPI_RECV(LL2fromSlaves(1),MAfSL_NPseudoExp,MPI_REAL,islave,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
        write(110,*)LL1fromSlaves(8)
        write(110,*) LL1fromSlaves(1:MAfSL_NPseudoExp)
!        write(110,*) LL2fromSlaves(1:MAfSL_NPseudoExp)
        
     enddo

  else        ! if one of the slaves
     write(120,*) 'receiving from master'
     call MPI_RECV(SL_InitPseudoExp,1,MPI_INT,root_process,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)      
     call MPI_RECV(SL_NPseudoExp,1,MPI_INT,root_process,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)      
     SL_FinalPseudoExp=SL_InitPseudoExp+SL_NPseudoExp-1
     write(120,*) SL_InitPseudoExp,SL_FinalPseudoExp,SL_NPseudoExp

     call random_seed()

! -----------------------------------------------------------------
! 3.1 Generate Poisson distribution about expected null value in each bin
! -----------------------------------------------------------------
     PlotObsEvts=0
     do iHypothesis=1,2
        !     do iPseudoExp=1,NPseudoExp
        do iPseudoExp=SL_InitPseudoExp,SL_FinalPseudoExp
!!RRremove1221           !        if( mod(iPseudoExp,10000).eq.0 ) print *, "Pseudo experiment ",iPseudoExp,"/",NPseudoExp
!!RRremove1221           GotNumEvents=.false.
!!RRremove1221           do iBin=1,NumBins
!!RRremove1221              do while (.not. GotNumEvents(iBin))
!!RRremove1221                 call random_number(nran(1:2))
!!RRremove1221                 nran(1)=nran(1)*PoissonMax(iHypothesis,iBin)        
!!RRremove1221                 TryEvts=int(5d0*ExpectedEvents(iHypothesis,iBin)*nran(2))
!!RRremove1221                 if (Poisson(ExpectedEvents(iHypothesis,iBin),TryEvts) .gt. nran(1)) then
!!RRremove1221                    ObsEvents(iHypothesis,iBin)=TryEvts
!!RRremove1221                    GotNumEvents(iBin)=.true.
!!RRremove1221                 endif
!!RRremove1221                 if (ExpectedEvents(iHypothesis,iBin) .eq. 0) then
!!RRremove1221                    ObsEvents(iHypothesis,iBin)=0
!!RRremove1221                    GotNumEvents(iBin)=.true.
!!RRremove1221                 endif
!!RRremove1221              enddo
!!RRremove1221              if (SUA .eq. 1) then
!!RRremove1221                 call random_number(sran)
!!RRremove1221                 sran=sran*(DeltaN**2-1d0)/DeltaN + 1d0/DeltaN
!!RRremove1221                 ObsEvents(iHypothesis,iBin)=ObsEvents(iHypothesis,iBin)*sran
!!RRremove1221              endif
!!RRremove1221              
!!RRremove1221              if (iBin .eq. 1) then
!!RRremove1221                 PlotObsEvts(iHypothesis,ObsEvents(iHypothesis,iBin))=PlotObsEvts(iHypothesis,ObsEvents(iHypothesis,iBin))+1
!!RRremove1221              endif
!!RRremove1221              
!!RRremove1221           enddo
!!RRremove1221           
!!RRremove1221           ! -----------------------------------------------------------------! 3.1a  change number of observed events in all bins for scale uncertainty    
!!RRremove1221           ! -----------------------------------------------------------------
!!RRremove1221           if (SUA .eq. 2) then
!!RRremove1221              call random_number(sran)
!!RRremove1221              sran=sran*(DeltaN**2-1d0)/DeltaN + 1d0/DeltaN
!!RRremove1221              write(201,*) sran
!!RRremove1221              ObsEvents(iHypothesis,1:NumBins)=ObsEvents(iHypothesis,1:NumBins)*sran
!!RRremove1221           endif
!!RRremove1221           ! -----------------------------------------------------------------
!!RRremove1221           ! 3.2 Now find the log likelihood, with a Poisson distr in each bin
!!RRremove1221           ! -----------------------------------------------------------------
!!RRremove1221           
!!RRremove1221           LLRatio=0d0
!!RRremove1221           do iBin=1,NumBins
!!RRremove1221              if (ObsEvents(iHypothesis,iBin) .ne. 0 .and. ExpectedEvents(1,iBin) .ne. 0 .and. &
!!RRremove1221                   & ExpectedEvents(2,iBin) .ne. 0) then
!!RRremove1221                 LLRatio=LLRatio+ObsEvents(iHypothesis,iBin)*dlog(1d0*ExpectedEvents(1,iBin)/ExpectedEvents(2,iBin))
!!RRremove1221              endif
!!RRremove1221           enddo
!!RRremove1221           
!!RRremove1221           ! offset to get the distributions with the histogram limits
!!RRremove1221           if (iHypothesis .eq. 1 .and. iPseudoExp .eq. 1) then
!!RRremove1221              offset=-100d0*(int(LLRatio)/100)
!!RRremove1221              print *, 'using offset of ', offset
!!RRremove1221           endif
!!RRremove1221           
!!RRremove1221           LLRatio=LLRatio+offset
!!RRremove1221           
!!RRremove1221           if (iPseudoExp .eq. 1) then
!!RRremove1221              print *, 'LLRatio=',LLRatio
!!RRremove1221           endif
!!RRremove1221
!! RR remove -- MPI debug
           LLRatio=1d0*iPseudoExp
           SL_LLRatio_array(iHypothesis,iPseudoExp)=LLRatio
           
!           print *, 'sending to master'
!           print *, iPseudoExp, SlavesLLRatio_array(1:2, iPseudoExp)
!           call MPI_SEND(iPseudoExp,1,MPI_INT,root_process,iPseudotag,MPI_COMM_WORLD,ierr)
!           call MPI_SEND(SlavesLLRatio_array(1:2,iPseudoExp),2,MPI_REAL,root_process,LLtag,MPI_COMM_WORLD,ierr)
        enddo

        
        do i=1,1000
           s=101+iHypothesis
           write(s,*) i,PlotObsEvts(iHypothesis,i)
        enddo
     enddo
     SL_LL1Ratio_array(SL_InitPseudoExp:SL_FinalPseudoExp)=SL_LLRatio_array(1,SL_InitPseudoExp:SL_FinalPseudoExp)
     SL_LL2Ratio_array(SL_InitPseudoExp:SL_FinalPseudoExp)=SL_LLRatio_array(1,SL_InitPseudoExp:SL_FinalPseudoExp)
     SL_inittag=1001
     SL_finaltag=1002
     LLtag1=10003
     LLtag2=10004
     write(120,*), 'sending to master from slave'
     write(120,*), SL_InitPseudoExp, SL_FinalPseudoExp
!     call MPI_SEND(SL_InitPseudoExp,1,MPI_INT,root_process,SL_inittag,MPI_COMM_WORLD,ierr)
!     call MPI_SEND(SL_NPseudoExp,1,MPI_INT,root_process,SL_inittag,MPI_COMM_WORLD,ierr)
!     call MPI_SEND(SL_LL1Ratio_array(SL_InitPseudoExp),SL_NPseudoExp,MPI_REAL, &
!          &root_process,LLtag1,MPI_COMM_WORLD,ierr)
!     call MPI_SEND(SL_LL2Ratio_array(SL_InitPseudoExp),SL_NPseudoExp,MPI_REAL, &
!          &root_process,LLtag2,MPI_COMM_WORLD,ierr)
     call MPI_SEND(SL_LL1Ratio_array(SL_InitPseudoExp),SL_NPseudoExp,MPI_REAL, &
          &root_process,SL_inittag,MPI_COMM_WORLD,ierr)
!    call MPI_SEND(SL_LL2Ratio_array(SL_InitPseudoExp),SL_NPseudoExp,MPI_REAL, &
!         &root_process,SL_inittag,MPI_COMM_WORLD,ierr)
     write(120,*), SL_LL1Ratio_array(SL_InitPseudoExp:SL_FinalPseudoExp)
     write(120,*), SL_LL2Ratio_array(SL_InitPseudoExp:SL_FinalPseudoExp)
     write(120,*), LL1_int(SL_InitPseudoExp:SL_FinalPseudoExp)
     

  endif

  
!!RRremove1221  if (my_id .eq.root_process) then
!!RRremove1221
!!RRremove1221     LLRatio_max=-1d-6
!!RRremove1221     LLRatio_min=1d6
!!RRremove1221     
!!RRremove1221     do iHypothesis=1,2
!!RRremove1221        do iPseudoExp=1,NPseudoExp
!!RRremove1221           if (LLRatio_array(iHypothesis,iPseudoExp) .gt. LLRatio_max) then
!!RRremove1221              LLRatio_max=LLRatio_array(iHypothesis,iPseudoExp) 
!!RRremove1221           endif
!!RRremove1221           if (LLRatio_array(iHypothesis,iPseudoExp) .lt. LLRatio_min) then
!!RRremove1221              LLRatio_min=LLRatio_array(iHypothesis,iPseudoExp) 
!!RRremove1221           endif
!!RRremove1221        enddo
!!RRremove1221     enddo
!!RRremove1221     
!!RRremove1221     LLRatio_max=ceiling(LLRatio_max)
!!RRremove1221     LLRatio_min=floor(LLRatio_min)
!!RRremove1221     
!!RRremove1221     if (LLRatio_max .lt. LLRatio_min) then
!!RRremove1221        print *, 'ERROR: MAX OF LL RATIO IS SMALLER THAN MIN!'
!!RRremove1221        print *, LLRatio_max,LLRatio_min
!!RRremove1221        stop
!!RRremove1221     endif
!!RRremove1221     if (LLRatio_min .eq. LLRatio_max) then
!!RRremove1221        LLRatio_max=LLRatio_max+1d0
!!RRremove1221        LLRatio_min=LLRatio_min-1d0
!!RRremove1221     endif
!!RRremove1221     print *, LLRatio_min,LLRatio_max
!!RRremove1221     
!!RRremove1221     LLHisto(1)%NBins   = 1000
!!RRremove1221     LLHisto(1)%BinSize = (LLRatio_max-LLRatio_min)/LLHisto(1)%NBins
!!RRremove1221     LLHisto(1)%LowVal  = LLRatio_min
!!RRremove1221     LLHisto(1)%Hits(:) = 0
!!RRremove1221     print *, 'using ', LLHisto(1)%BinSize, ' bins in range ',LLRatio_min, LLRatio_max
!!RRremove1221     LLHisto(2)%NBins   = LLHisto(1)%NBins
!!RRremove1221     LLHisto(2)%BinSize = LLHisto(1)%BinSize
!!RRremove1221     LLHisto(2)%LowVal  = LLHisto(1)%LowVal
!!RRremove1221     LLHisto(2)%Hits(:) = LLHisto(1)%Hits(:)
!!RRremove1221     
!!RRremove1221     do iHypothesis=1,2
!!RRremove1221        do iPseudoExp=1,NPseudoExp
!!RRremove1221           LLRatio=LLRatio_array(iHypothesis,iPseudoExp)
!!RRremove1221           WhichBin = (LLRatio-LLHisto(iHypothesis)%LowVal)/LLHisto(iHypothesis)%BinSize + 1
!!RRremove1221           WhichBin=int(WhichBin)
!!RRremove1221           if( WhichBin.lt.0 ) WhichBin = 1
!!RRremove1221           if( WhichBin.gt.LLHisto(iHypothesis)%NBins ) WhichBin = LLHisto(iHypothesis)%NBins
!!RRremove1221           LLHisto(iHypothesis)%Hits(WhichBin) = LLHisto(iHypothesis)%Hits(WhichBin) + 1
!!RRremove1221        enddo
!!RRremove1221        
!!RRremove1221        ! -----------------------------------------------------------------
!!RRremove1221        ! 4. Calculate the integral of the log likelihood ratio curve, at each bin
!!RRremove1221        ! -----------------------------------------------------------------
!!RRremove1221        
!!RRremove1221        
!!RRremove1221        do LLbin=1,LLHisto(iHypothesis)%NBins
!!RRremove1221           IntLLRatio(iHypothesis)=IntLLRatio(iHypothesis)+&
!!RRremove1221                LLHisto(iHypothesis)%BinSize * LLHisto(iHypothesis)%Hits(LLbin)
!!RRremove1221           alpha(iHypothesis,LLBin)=IntLLRatio(iHypothesis)/(NPseudoExp*LLHisto(iHypothesis)%BinSize)
!!RRremove1221        enddo
!!RRremove1221     enddo
!!RRremove1221     print *, ""
!!RRremove1221     print *, "Writing log-likelihood distribution to output file"
!!RRremove1221     print *, ""
!!RRremove1221     do LLbin=1,LLHisto(1)%NBins
!!RRremove1221        !         print *, LLBin
!!RRremove1221        write(14,"(2X,I4,2X,1PE16.8,2X,I10,2X,1PE16.8,2X,I10,2X,1PE16.8,2X,1PE16.8)") LLbin, LLHisto(1)%LowVal+LLbin*LLHisto(1)%BinSize, LLHisto(1)%Hits(LLbin),LLHisto(2)%LowVal+LLbin*LLHisto(2)%BinSize, LLHisto(2)%Hits(LLbin),alpha(1,LLBin),alpha(2,LLBin)
!!RRremove1221     enddo
!!RRremove1221     
!!RRremove1221     
!!RRremove1221     !************************************************************
!!RRremove1221     !  5. Now find the point at which alpha=1-beta
!!RRremove1221     !************************************************************
!!RRremove1221     
!!RRremove1221     ! for the time being, I'm going to assume that both distributions have the same range and binning - can change this later
!!RRremove1221     do LLbin=1,LLHisto(1)%NBins
!!RRremove1221        checkalpha(LLBin)=alpha(1,LLBin)+alpha(2,LLBin)-1d0
!!RRremove1221        !   print *, LLBIn,alpha(1,LLBin),alpha(2,LLBin),checkalpha(LLBin)
!!RRremove1221        if (checkalpha(LLBin)*checkalpha(LLBin-1) .lt. 0d0) then
!!RRremove1221           alphamin=alpha(1,LLBin-1)
!!RRremove1221           alphamax=alpha(1,LLBin)
!!RRremove1221           betamin =alpha(2,LLBin-1)
!!RRremove1221           betamax =alpha(2,LLBin)
!!RRremove1221           
!!RRremove1221           if (checkalpha(LLBin) .eq. 1d0 .and. checkalpha(LLBin-1) .eq.-1d0) then
!!RRremove1221              ! this is comparing two identical hypotheses!
!!RRremove1221              alphamin=0.5d0
!!RRremove1221              alphamax=0.5d0
!!RRremove1221              betamin =0.5d0
!!RRremove1221              betamax =0.5d0
!!RRremove1221           endif
!!RRremove1221        endif
!!RRremove1221     enddo
!!RRremove1221
!!RRremove1221
!!RRremove1221     alpha_ave=(alphamin+alphamax)/2d0
!!RRremove1221     beta_ave=(betamin+betamax)/2d0
!!RRremove1221     sigs=dsqrt(2d0)*inverf(1d0-alpha_ave)
!!RRremove1221     
!!RRremove1221     print *, 'alpha value in range:', alphamin,alphamax
!!RRremove1221     print *, 'beta value in range:', betamin,betamax 
!!RRremove1221     print *, 'alpha average :', alpha_ave
!!RRremove1221     print *, 'beta average :', beta_ave
!!RRremove1221     print *, 'sigma value :', sigs
!!RRremove1221     
!!RRremove1221     write(14,"(A,2X,1PE16.8,2X,1PE16.8)") "# alpha value in range:",alphamin,alphamax
!!RRremove1221     write(14,"(A,2X,1PE16.8,2X,1PE16.8)") "# beta value in range:",betamin,betamax
!!RRremove1221     write(14,"(A,2X,1PE16.8,2X,1PE16.8)") "# alpha average:", alpha_ave 
!!RRremove1221     write(14,"(A,2X,1PE16.8,2X,1PE16.8)") "# beta average:", beta_ave   
!!RRremove1221     write(14,"(A,2X,1PE16.8,2X,1PE16.8)") "# sigma value: ", sigs
!!RRremove1221     
!!RRremove1221     
!!RRremove1221     write(15,"(1PE16.8)") alpha_ave
!!RRremove1221     write(15,"(1PE16.8)") beta_ave
!!RRremove1221     write(15,"(1PE16.8)") sigs
!!RRremove1221  endif

  call MPI_FINALIZE(ierr)

  close(12)
  close(13)
  close(112)
  close(113)
  close(14)
  close(15)
  stop


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


end program BinnedLogL
  
