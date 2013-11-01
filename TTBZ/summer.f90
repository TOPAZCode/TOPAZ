PROGRAM summer
use ifport
implicit none
integer :: NumArgs,NArg,op,iHisto,nrebin,NumEvents,NumPseudoExp,Lumi
character :: filename_in*(80),filename_in2*(80),filename_out*(50),operation*(10),factor_str*(20)
character :: iHisto_str*(5),nrebin_str*(5),NumEvents_str*(6),Lumi_str*(6),NumPseudoExp_str*(9)
real(8):: factor =1d10

! get number of arguments
  NumArgs = NArgs()-1

! get mode of operation
  call GetArg(1,operation)
  if( trim(operation).eq.'add' ) then
      op=1
  elseif( trim(operation).eq.'avg' ) then
      op=2
  elseif( trim(operation).eq.'sumbins' ) then
      op=3
  elseif( trim(operation).eq.'mul' ) then
      op=4
  elseif( trim(operation).eq.'quo' ) then
      op=5
  elseif( trim(operation).eq.'rebin' ) then
      op=6
  elseif( trim(operation).eq.'readlhe' ) then
      op=7
  elseif( trim(operation).eq.'likelihood' ) then
      op=8
  else
      write(*,*) "operation not available (add,avg,sumbins,mul,quo,rebin,readlhe,likelihood)",trim(operation)
      stop
  endif

! check for invalid arguments
  if( (op.lt.3 .and. NumArgs.le.3) .or. (op.eq.3 .and. NumArgs.ne.3) .or. (op.eq.4 .and. NumArgs.ne.4) .or. (op.eq.5 .and. NumArgs.le.3)  & 
                                   .or. (op.eq.6 .and. NumArgs.ne.5) .or. (op.eq.7 .and. NumArgs.ne.3) .or. (op.eq.8 .and. NumArgs.ne.7)  &
                                   .or. NumArgs.gt.49 ) then
    print *, "invalid number of arguments",op,NumArgs
    stop
  endif


  if( op.eq.1 ) then
      call GetArg(NumArgs,filename_out)
      open(unit=10,file=trim(filename_out),form='formatted',access='sequential',status='replace') ! open output file (last one in arguments list)
      print *, "adding histograms"
      do NArg=2,NumArgs-1  ! open input files
        call GetArg(NArg,filename_in)
        open(unit=10+NArg,file=trim(filename_in),form='formatted',access='sequential')
      enddo
      call add_avg(NumArgs,op)
      print *, "output written to "//trim(filename_out)

  elseif( op.eq.2 ) then
      call GetArg(NumArgs,filename_out)
      open(unit=10,file=trim(filename_out),form='formatted',access='sequential',status='replace') ! open output file (last one in arguments list)
      print *, "calculating mean value of histograms"
      do NArg=2,NumArgs-1  ! open input files
        call GetArg(NArg,filename_in)
        open(unit=10+NArg,file=trim(filename_in),form='formatted',access='sequential')
      enddo
      call add_avg(NumArgs,op)
      print *, "output written to "//trim(filename_out)

  elseif( op.eq.3 ) then
      print *, "summing up histogram bins"
      call GetArg(2,filename_in)
      open(unit=12,file=trim(filename_in),form='formatted',access='sequential')  ! open input file
      call GetArg(3,iHisto_str)
      read(iHisto_str,"(I2)") iHisto
      call sum_histobins(NumArgs,iHisto)

  elseif( op.eq.4 ) then
      call GetArg(NumArgs,filename_out)
      open(unit=10,file=trim(filename_out),form='formatted',access='sequential',status='replace') ! open output file (last one in arguments list)
      call GetArg(3,filename_in)
      open(unit=12,file=trim(filename_in),form='formatted',access='sequential')  ! open input file
      call GetArg(2,factor_str)
      read(factor_str,"(F8.6)") factor
      print *, "multiply histogram with ",factor
      call multiply_histo(NumArgs,factor)
      print *, "output written to "//trim(filename_out)

  elseif( op.eq.5 ) then
      call GetArg(NumArgs,filename_out)
      open(unit=10,file=trim(filename_out),form='formatted',access='sequential',status='replace') ! open output file (last one in arguments list)
      print *, "dividing histograms"
      do NArg=2,NumArgs-1  ! open input files
        call GetArg(NArg,filename_in)
        open(unit=10+NArg,file=trim(filename_in),form='formatted',access='sequential')
      enddo
      call divide(NumArgs,op)
      print *, "output written to "//trim(filename_out)

  elseif( op.eq.6 ) then
      print *, "rebinning"
      call GetArg(2,filename_in)
      open(unit=12,file=trim(filename_in),form='formatted',access='sequential')  ! open input file
      call GetArg(3,filename_out)
      open(unit=13,file=trim(filename_out),form='formatted',access='sequential')  ! open output file
      call GetArg(4,iHisto_str)
      read(iHisto_str,"(I2)") iHisto
      call GetArg(5,nrebin_str)
      read(nrebin_str,"(I2)") nrebin
      call rebinning(iHisto,nrebin)

  elseif( op.eq.7 ) then
      print *, "reading lhe event file"
      call GetArg(2,filename_in)
      open(unit=12,file=trim(filename_in),form='formatted',access='sequential')  ! open input file
      call GetArg(3,filename_out)
      open(unit=13,file=trim(filename_out),form='formatted',access='sequential')  ! open output file
      call readlhe()

  elseif( op.eq.8 ) then! SYNTAX: summer likelihood infile1 infile2 Histo Events PseudoEvents outfile
      call GetArg(2,filename_in)
      open(unit=12,file=trim(filename_in),form='formatted',access='sequential')  ! open input file1
      call GetArg(3,filename_in2)
      open(unit=13,file=trim(filename_in2),form='formatted',access='sequential')  ! open input file2
      call GetArg(4,iHisto_str)
      read(iHisto_str,"(I2)") iHisto   

!       call GetArg(5,NumEvents_str)
!       read(NumEvents_str,"(I6)") NumEvents
      call GetArg(5,Lumi_str)
      read(Lumi_str,"(I6)") Lumi

      call GetArg(6,NumPseudoExp_str)
      read(NumPseudoExp_str,"(I9)") NumPseudoExp
      call GetArg(7,filename_out)
      open(unit=14,file=trim(filename_out),form='formatted',access='sequential')  ! open output file
      write(*,*) ""
!       write(*,"(A,I2,A,I6,A,I9,A)") "reading histogram ",iHisto," for likelihood analysis with ",NumEvents," events using ",NumPseudoExp," pseudo-experiments"
      write(*,"(A,I2,A,I6,A,I9,A)") "reading histogram ",iHisto," for likelihood analysis with Lumi=",Lumi,"fb^-1 using ",NumPseudoExp," pseudo-experiments"
      write(*,*) ""
!       call likelihood(iHisto,NumEvents,NumPseudoExp)
      call likelihood(iHisto,Lumi,NumPseudoExp)
  endif


! close all files
  do NArg=2,NumArgs
      close(unit=10+NArg)
  enddo
  close(10)

END PROGRAM






SUBROUTINE add_avg(NumArgs,op)
implicit none
integer :: NArg,row,NumArgs,op
character(len=*),parameter :: fmt1 = "(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)"
integer,parameter :: MaxFiles=50
integer :: NHisto(1:MaxFiles)=-999999,Hits(1:MaxFiles)=-999999
real(8) :: BinVal(1:MaxFiles)=-1d-99,Value(1:MaxFiles)=-1d-99,Error(1:MaxFiles)=-1d-99
real(8) :: SumValue,SumError,n
integer :: SumHits
character :: dummy*(1)


  do while(.not.eof(12))  ! loop over all rows
    SumValue = 0d0
    SumError = 0d0
    Sumhits  = 0

    do NArg=1,NumArgs-2  ! loop over all input files
      read(unit=11+NArg,fmt="(A)") dummy
      if(dummy(1:1).eq."#") cycle
      backspace(unit=11+NArg) ! go to the beginning of the line
      read(unit=11+NArg,fmt=fmt1) NHisto(NArg),dummy,BinVal(NArg),dummy,Value(NArg),dummy,Error(NArg),dummy,Hits(NArg),dummy
      n = dble(NumArgs-2)
      if( op.eq.1 ) then
          SumValue = SumValue + Value(NArg)
          SumError = SumError + Error(NArg)**2
          Sumhits  = Sumhits + Hits(NArg)
      elseif( op.eq.2 ) then
          SumValue = SumValue + Value(NArg) / n
          Sumhits  = Sumhits + Hits(NArg)
      endif
    enddo

    if( op.eq.2 ) then
        do NArg=1,NumArgs-2
            SumError = SumError + (SumValue-Value(NArg))**2 /n /(n-1d0)
        enddo
    endif

    if(dummy(1:1).eq."#") cycle
    if(.not. all(NHisto(1:NumArgs-2).eq.NHisto(1))) print *, "Error: unequal values for NHisto in row ",row
    if(.not. all(BinVal(1:NumArgs-2).eq.BinVal(1))) print *, "Error: unequal values for BinVal in row ",row

    SumError = dsqrt(SumError)

    if(abs(SumHits).gt.999999999) SumHits=999999999
    write(* ,fmt1) NHisto(1),"|",BinVal(1),"|",SumValue,"|",SumError,"|",SumHits,"|"
    write(10,fmt1) NHisto(1),"|",BinVal(1),"|",SumValue,"|",SumError,"|",SumHits,"|"
  enddo

END SUBROUTINE




SUBROUTINE divide(NumArgs,op)
implicit none
integer :: NArg,row,NumArgs,op
character(len=*),parameter :: fmt1 = "(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)"
integer,parameter :: MaxFiles=50
integer :: NHisto(1:MaxFiles)=-999999,Hits(1:MaxFiles)=-999999
real(8) :: BinVal(1:MaxFiles)=-1d-99,Value(1:MaxFiles)=-1d-99,Error(1:MaxFiles)=-1d-99
real(8) :: SumValue,SumError,n
integer :: SumHits
character :: dummy*(1)


  do while(.not.eof(12))  ! loop over all rows
    SumValue = 0d0
    SumError = 0d0
    Sumhits  = 0

    do NArg=1,NumArgs-2  ! loop over all input files
      read(unit=11+NArg,fmt="(A)") dummy
      if(dummy(1:1).eq."#") cycle
      backspace(unit=11+NArg) ! go to the beginning of the line
      read(unit=11+NArg,fmt=fmt1) NHisto(NArg),dummy,BinVal(NArg),dummy,Value(NArg),dummy,Error(NArg),dummy,Hits(NArg),dummy
      n = dble(NumArgs-2)
    enddo


    if(dummy(1:1).eq."#") cycle
    if(.not. all(NHisto(1:NumArgs-2).eq.NHisto(1))) print *, "Error: unequal values for NHisto in row ",row
    if(.not. all(BinVal(1:NumArgs-2).eq.BinVal(1))) print *, "Error: unequal values for BinVal in row ",row

    SumError = dsqrt((Error(1)/Value(2))**2+(Error(2)*Value(1))**2)
    if(SumHits.gt.999999999) SumHits=999999999
    write(* ,fmt1) NHisto(1),"|",BinVal(1),"|",Value(1)/Value(2),"|",SumError,"|",SumHits,"|"
    write(10,fmt1) NHisto(1),"|",BinVal(1),"|",Value(1)/Value(2),"|",SumError,"|",SumHits,"|"
  enddo

END SUBROUTINE




SUBROUTINE sum_histobins(NumArgs,iHisto)
implicit none
character(len=*),parameter :: fmt1 = "(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)"
integer :: NumArgs,NArg,iHisto
integer,parameter :: MaxFiles=50
integer :: NHisto(1:MaxFiles)=-999999,Hits(1:MaxFiles)=-999999
real(8) :: BinVal(1:MaxFiles)=-1d-99,Value(1:MaxFiles)=-1d-99,Error(1:MaxFiles)=-1d-99
real(8) :: SumValue,BinSize_Tmp,BinSize
character :: dummy*(1)


  BinSize_Tmp = 1d-100
  BinSize = 1d-100

  SumValue = 0d0
  do while(.not.eof(12))  ! loop over all rows

      NArg=1
      read(unit=11+NArg,fmt="(A)") dummy
      if(dummy(1:1).eq."#") cycle
      backspace(unit=11+NArg) ! go to the beginning of the line
      read(unit=11+NArg,fmt=fmt1) NHisto(1),dummy,BinVal(1),dummy,Value(1),dummy,Error(1),dummy,Hits(1),dummy

      if(NHisto(1).ne.iHisto) cycle

      ! detect bin size (assume same size for all bins)
      if(BinSize_Tmp.ne.1d-100 .and. BinSize.eq.1d-100) then
          BinSize=dabs(BinVal(1)-BinSize_Tmp)
!           print *, "Bin size=",BinSize
      endif
      if(BinSize_Tmp.eq.1d-100) BinSize_Tmp=BinVal(1)

      SumValue = SumValue  + Value(1)
!       write(* ,fmt1) NHisto(1),"|",BinVal(1),"|",Value(1),"|",Error(1),"|",Hits(1),"|"

  enddo

  SumValue = SumValue * BinSize
  write(*,"(1PE23.16)") SumValue
!   write(10,"(1PE23.16)") SumValue

END SUBROUTINE


SUBROUTINE multiply_histo(NumArgs,factor)
implicit none
character(len=*),parameter :: fmt1 = "(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)"
integer :: NumArgs,NArg
integer,parameter :: MaxFiles=50
integer :: NHisto(1:MaxFiles)=-999999,Hits(1:MaxFiles)=-999999
real(8) :: BinVal(1:MaxFiles)=-1d-99,Value(1:MaxFiles)=-1d-99,Error(1:MaxFiles)=-1d-99
real(8) :: SumValue,factor
character :: dummy*(1)

  do while(.not.eof(12))  ! loop over all rows
      NArg=1
      read(unit=11+NArg,fmt="(A)") dummy
      if(dummy(1:1).eq."#") cycle
      backspace(unit=11+NArg) ! go to the beginning of the line
      read(unit=11+NArg,fmt=fmt1) NHisto(NArg),dummy,BinVal(NArg),dummy,Value(NArg),dummy,Error(NArg),dummy,Hits(NArg),dummy

      if(dummy(1:1).eq."#") cycle
      if(Hits(1).gt.999999999) Hits=999999999
      write(* ,fmt1) NHisto(1),"|",BinVal(1),"|",factor*Value(1),"|",factor*Error(1),"|",Hits(1),"|"
      write(10,fmt1) NHisto(1),"|",BinVal(1),"|",factor*Value(1),"|",factor*Error(1),"|",Hits(1),"|"
  enddo


END SUBROUTINE



SUBROUTINE rebinning(iHisto,nrebin)
implicit none
character(len=*),parameter :: fmt1 = "(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)"
integer :: NumArgs,NArg,iHisto,nrebin, ibin
integer,parameter :: MaxFiles=50
integer :: NHisto(1:MaxFiles)=-999999,Hits(1:MaxFiles)=-999999
real(8) :: BinVal(1:MaxFiles)=-1d-99,Value(1:MaxFiles)=-1d-99,Error(1:MaxFiles)=-1d-99
real(8) :: newvalue, newerror
integer :: newhits
real(8) :: SumValue,BinSize_Tmp,BinSize
character :: dummy*(1)


  BinSize_Tmp = 1d-100
  BinSize = 1d-100

  SumValue = 0d0
  do while(.not.eof(12))  ! loop over all rows


      read(unit=12,fmt="(A)") dummy
      if(dummy(1:1).eq."#") cycle
      backspace(unit=12) ! go to the beginning of the line
      read(unit=12,fmt=fmt1) NHisto(1),dummy,BinVal(1),dummy,Value(1),dummy,Error(1),dummy,Hits(1),dummy

      if(NHisto(1).eq.iHisto) then
!         write(* ,fmt1) NHisto(1),"|",BinVal(1),"|",Value(1),"|",Error(1),"|",Hits(1),"|"

         read(unit=12,fmt="(A)") dummy
         if(dummy(1:1).eq."#") cycle
         backspace(unit=12) ! go to the beginning of the line

         newvalue = Value(1) 
         newerror = Error(1)**2
         newhits = Hits(1)
         
         do ibin = 1, nrebin
            read(unit=12,fmt=fmt1) NHisto(1+ibin),dummy,BinVal(1+ibin),dummy,Value(1+ibin),dummy,Error(1+ibin),dummy,Hits(1+ibin),dummy
            !         write(* ,fmt1) NHisto(2),"|",BinVal(2),"|",Value(2),"|",Error(2),"|",Hits(2),"|"
            newvalue = newvalue + Value(1+ibin)
            newerror = newerror + Error(1+ibin)**2
            newhits = newhits + Hits(1+ibin)
         enddo

         newerror = dsqrt(newerror)/real(nrebin+1)
         newvalue = newvalue/real(nrebin+1)

         write(13 ,fmt1) NHisto(1),"|",BinVal(1),"|",newvalue,"|",newerror,"|",newhits,"|"
         do ibin = 1, nrebin
            write(13 ,fmt1) NHisto(1),"|",BinVal(1+ibin),"|",newvalue,"|",newerror,"|",newhits,"|"
         enddo

         write(*  ,fmt1) NHisto(1),"|",BinVal(1),"|",newvalue,"|",newerror,"|",newhits,"|"

      else
         write(13 ,fmt1) NHisto(1),"|",BinVal(1),"|",Value(1),"|",Error(1),"|",Hits(1),"|"
         write(*  ,fmt1) NHisto(1),"|",BinVal(1),"|",Value(1),"|",Error(1),"|",Hits(1),"|"
      endif

  enddo


END SUBROUTINE




SUBROUTINE readlhe()
implicit none
character(len=*),parameter :: fmt1 = "(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)"
integer :: NumArgs,NArg
integer,parameter :: MaxFiles=50
integer :: NHisto(1:MaxFiles)=-999999,Hits(1:MaxFiles)=-999999
real(8) :: BinVal(1:MaxFiles)=-1d-99,Value(1:MaxFiles)=-1d-99,Error(1:MaxFiles)=-1d-99
real(8) :: SumValue,factor
character :: dummy*(1)


  do while(.not.eof(12))  ! loop over all events

      read(unit=12,fmt="(A)") dummy
      if(dummy(1:1).eq."#") cycle
      backspace(unit=12) ! go to the beginning of the line
      read(unit=12,fmt=fmt1) NHisto(NArg),dummy,BinVal(NArg),dummy,Value(NArg),dummy,Error(NArg),dummy,Hits(NArg),dummy

  enddo

END SUBROUTINE





SUBROUTINE likelihood(iHisto,LumiORNumEvents,NumPseudoExp)
implicit none
integer :: iHisto,LumiORNumEvents,NumPseudoExp
integer,parameter :: MaxBins=1000, MaxEvents=200000
integer,parameter :: NumLLbins=1000! the LLratio should only vary between - and + this number
character(len=*),parameter :: fmt1 = "(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)"
character :: dummy*(1)
integer :: iBin,iPseudoExp,NumBins1=0,NumBins2=0,ranBin,BinAccept(1:MaxEvents),NumAccepted,NEvent
integer :: NumEvents,NumExpectedEvents(1:2),iHypothesis
integer :: NHisto(1:2)=-999999,Hits(1:2,1:MaxBins)=-999999
real(8) :: Lumi,BinVal(1:2,1:MaxBins)=-1d-99,Value(1:2,1:MaxBins)=-1d-99,Error(1:2,1:MaxBins)=-1d-99
real(8) :: ymax,zran(1:2),ValAccepted(1:MaxEvents),LLRatio,WhichBin,BinSize
real(8) :: sigmatot(1:2),check(1:2)
integer :: SelectedEvent,LLbin
type :: Histogram
    integer :: NBins
    real(8) :: BinSize
    real(8) :: LowVal
    integer :: Hits(1:500)
end type
type(Histogram) :: LLHisto(1:2)
logical, parameter :: normalize=.true.



! NumEvents=LumiORNumEvents
Lumi = dble(LumiORNumEvents)


!--------------------------------------------------
!          0. init
!--------------------------------------------------
LLHisto(1)%NBins   = 1000
LLHisto(1)%BinSize = 1.0d0
LLHisto(1)%LowVal  = -500d0
LLHisto(1)%Hits(:) = 0

LLHisto(2)%NBins   = LLHisto(1)%NBins
LLHisto(2)%BinSize = LLHisto(1)%BinSize
LLHisto(2)%LowVal  = LLHisto(1)%LowVal
LLHisto(2)%Hits(:) = LLHisto(1)%Hits(:)



!--------------------------------------------------
!          1. reading input files
!--------------------------------------------------
  do while(.not.eof(12))!   reading input file 1
      read(unit=12,fmt="(A)") dummy
      if(dummy(1:1).eq."#") cycle
      backspace(unit=12) ! go to the beginning of the line

      read(unit=12,fmt=fmt1) NHisto(1),dummy,BinVal(1,NumBins1+1),dummy,Value(1,NumBins1+1),dummy,Error(1,NumBins1+1),dummy,Hits(1,NumBins1+1),dummy
      if( NHisto(1).ne.iHisto ) cycle
      NumBins1=NumBins1 + 1
  enddo

  do while(.not.eof(13))!   reading input file 2
      read(unit=13,fmt="(A)") dummy
      if(dummy(1:1).eq."#") cycle
      backspace(unit=13) ! go to the beginning of the line

      read(unit=13,fmt=fmt1) NHisto(2),dummy,BinVal(2,NumBins2+1),dummy,Value(2,NumBins2+1),dummy,Error(2,NumBins2+1),dummy,Hits(2,NumBins2+1),dummy
      if( NHisto(2).ne.iHisto ) cycle
      NumBins2=NumBins2 + 1
  enddo

  if( NumBins1.ne.NumBins2 ) then
     print *, "Error: Number of bins in input file 1 and 2 are different: ",NumBins1,NumBins2
     stop
  endif

  sigmatot(1:2) = 0d0
  BinSize = BinVal(1,2)-BinVal(1,1)
  do iBin=1,NumBins1! calculate total cross section
       sigmatot(1) = sigmatot(1) + Value(1,iBin) * BinSize
       sigmatot(2) = sigmatot(2) + Value(2,iBin) * BinSize
  enddo
  NumExpectedEvents(1) = int(Lumi * sigmatot(1) *8 ) ! factor of 8 from lepton species
  NumExpectedEvents(2) = int(Lumi * sigmatot(2) *8 ) 


! normalize distributions
  if( normalize ) then
   print *, "Normalizing input files" 
   do iBin=1,NumBins1! normalize histogram
       Value(1,iBin) = Value(1,iBin)*BinSize/sigmatot(1)
       Value(2,iBin) = Value(2,iBin)*BinSize/sigmatot(2)
   enddo
  endif

! print the input histograms 
  write(*,"(2X,A,16X,A,11X,A,16X,A)") "NBin|","Input file 1","|","Input file 2"
  do iBin=1,NumBins1
    write(*,fmt="(2X,1I3,A,2X,1PE10.3,2X,1PE23.16,A,2X,1PE10.3,2X,1PE23.16)") iBin," | ",BinVal(1,iBin),Value(1,iBin)," | ",BinVal(2,iBin),Value(2,iBin)
    if( dabs(BinVal(1,iBin)-BinVal(2,iBin)).gt.1d-6 ) then
        print *, "Error: Different bin sizes in input files 1 and 2"
        stop
    endif
  enddo
  write(*,"(A,1PE16.8,A,I6)") "Total cross section of input file 1: ",sigmatot(1),"   <-->   Number of events: ",NumExpectedEvents(1)
  write(*,"(A,1PE16.8,A,I6)") "Total cross section of input file 2: ",sigmatot(2),"   <-->   Number of events: ",NumExpectedEvents(2)
  check(1:2) = 0d0
  do iBin=1,NumBins1! check
      check(1) = check(1) + Value(1,iBin)
      check(2) = check(2) + Value(2,iBin)
  enddo
  write(*,"(A,1PE16.8)") "Sum of bins in input file 1: ",check(1)
  write(*,"(A,1PE16.8)") "Sum of bins in input file 2: ",check(2)
  write(*,*) ""









!----------------------------------------------------------------------------------------
!          2. generate pseudo experiments for null hypothesis and alternative hypothesis
!----------------------------------------------------------------------------------------
do iHypothesis=1,2

!   NumEvents = NumExpectedEvents(1)              ! choose number of events fixed for each hypothesis
  NumEvents = NumExpectedEvents(iHypothesis)    ! choose number of events differently for each hypothesis

  write(*,"(A,I6,A)") "Generating ",NumEvents," events"
  if( NumEvents.gt.MaxEvents ) then
    print *, "Error: NumEvents is too large. Increase MaxEvents in SUBROUTINE likelihood."
    stop
  endif


  do iPseudoExp=1,NumPseudoExp
      if( mod(iPseudoExp,1000).eq.0 ) print *, "Pseudo experiment ",iPseudoExp,"/",NumPseudoExp

!************************************************************
!  2.1 generate event sample according to hypothesis
!************************************************************
      ! find max. of y-axis
      ymax=-1d13;
      do iBin=1,NumBins1
         if(  Value(iHypothesis,iBin).gt.ymax ) ymax = Value(iHypothesis,iBin)
      enddo

      NumAccepted=0
      call random_seed()     
      do while( NumAccepted .lt. NumEvents )!  loop until required number of events is reached
           call random_number(zran(1:2))         
           ranBin = 1 + int( zran(1)*NumBins1 ) ! randomly select a bin (=event)
           if( ymax*zran(2) .lt.  Value(iHypothesis,ranBin) ) then! accept/reject this event
               NumAccepted=NumAccepted+1
               BinAccept(NumAccepted) = ranBin
               ValAccepted(NumAccepted) = Value(iHypothesis,ranBin)
           endif
      enddo


!************************************************************
!  2.2 calculate log likelihood ratio 
!************************************************************
      LLRatio = 0d0
      do NEvent=1,NumEvents
         SelectedEvent = BinAccept(NEvent)
         if( ValAccepted(NEvent).ne.Value(iHypothesis,SelectedEvent) ) then! can be removed/simplified later 
             print *, "error in ll"
             stop
         endif
         LLRatio = LLRatio + 2d0*dlog( (Value(1,SelectedEvent))/(Value(2,SelectedEvent)) )
      enddo


!************************************************************
!  2.3 bin the likelihood value 
!************************************************************
      WhichBin = (LLRatio-LLHisto(iHypothesis)%LowVal)/LLHisto(iHypothesis)%BinSize + 1
      if( WhichBin.lt.0 ) WhichBin = 1
      if( WhichBin.gt.LLHisto(iHypothesis)%NBins ) WhichBin = LLHisto(iHypothesis)%NBins
      LLHisto(iHypothesis)%Hits(WhichBin) = LLHisto(iHypothesis)%Hits(WhichBin) + 1

  enddo! iPseudoExp


!************************************************************
!  3. write out the LL distribution
!************************************************************
      print *, ""
      print *, "Writing log-likelihood distribution to outful file"
      print *, ""
      do LLbin=1,LLHisto(iHypothesis)%NBins
           write(14,"(2X,I1,2X,I4,2X,1PE16.8,I10)") 1,LLbin, LLHisto(iHypothesis)%LowVal+LLbin*LLHisto(iHypothesis)%BinSize, LLHisto(iHypothesis)%Hits(LLbin)
      enddo

enddo! iHypothesis




END SUBROUTINE
























