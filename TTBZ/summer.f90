PROGRAM summer
use ifport
implicit none
integer :: NumArgs,NArg,op,iHisto,nrebin,NumEvents,NumPseudoExp,Lumi
character :: filename_in*(80),filename_in2*(80),filename_out*(50),operation*(10),factor_str*(20)
character :: iHisto_str*(5),nrebin_str*(5),NumEvents_str*(6),Lumi_str*(6),NumPseudoExp_str*(9),DeltaN_str*(5)
character :: dVcoupl_str*(10),dAcoupl_str*(10)
real(8) :: DeltaN
real(8):: factor =1d10,dVcoupl(1:6)=1d10,dAcoupl(1:6)=1d10,VcouplSM,AcouplSM
logical :: silent=.false.


! get number of arguments
  NumArgs = NArgs()-1

! get mode of operation
  call GetArg(1,operation)
  if( trim(operation).eq.'add' ) then
      op=1
  elseif( trim(operation).eq.'addsilent' ) then
      op=1
      silent=.true.
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
  elseif( trim(operation).eq.'ttbzcoupl' ) then
      op=9
  else
      write(*,*) "operation not available (add,avg,sumbins,mul,quo,rebin,readlhe,likelihood)",trim(operation)
      stop
  endif

! check for invalid arguments
  if( (op.lt.3 .and. NumArgs.le.3) .or. (op.eq.3 .and. NumArgs.ne.3) .or. (op.eq.4 .and. NumArgs.ne.4) .or. (op.eq.5 .and. NumArgs.le.3)  & 
                                   .or. (op.eq.6 .and. NumArgs.ne.5) .or. (op.eq.7 .and. NumArgs.ne.3) .or. (op.eq.8 .and. NumArgs.ne.8)  &
                                   .or. (op.eq.9 .and. NumArgs.ne.21).or. NumArgs.gt.49 ) then
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
      call add_avg(NumArgs,op,silent)
      print *, "output written to "//trim(filename_out)

  elseif( op.eq.2 ) then
      call GetArg(NumArgs,filename_out)
      open(unit=10,file=trim(filename_out),form='formatted',access='sequential',status='replace') ! open output file (last one in arguments list)
      print *, "calculating mean value of histograms"
      do NArg=2,NumArgs-1  ! open input files
        call GetArg(NArg,filename_in)
        open(unit=10+NArg,file=trim(filename_in),form='formatted',access='sequential')
      enddo
      call add_avg(NumArgs,op,silent)
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

  elseif( op.eq.8 ) then! SYNTAX: summer likelihood infile1 infile2 Histo Events PseudoEvents DeltaN outfile
      call GetArg(2,filename_in)
      open(unit=12,file=trim(filename_in),form='formatted',access='sequential')  ! open input file1
      call GetArg(3,filename_in2)
      open(unit=13,file=trim(filename_in2),form='formatted',access='sequential')  ! open input file2
      print *, 'Using input files:', filename_in, filename_in2
      call GetArg(4,iHisto_str)
      read(iHisto_str,"(I2)") iHisto   

!       call GetArg(5,NumEvents_str)
!       read(NumEvents_str,"(I6)") NumEvents
      call GetArg(5,Lumi_str)
      read(Lumi_str,"(I6)") Lumi

      call GetArg(6,NumPseudoExp_str)
      call GetArg(7,DeltaN_str)
      read(DeltaN_str,"(F9.6)") DeltaN
!      DeltaN=1.3d0
      read(NumPseudoExp_str,"(I9)") NumPseudoExp
      call GetArg(8,filename_out)
      open(unit=14,file=trim(filename_out),form='formatted',access='sequential')  ! open output file
      write(*,*) ""
!       write(*,"(A,I2,A,I6,A,I9,A)") "reading histogram ",iHisto," for likelihood analysis with ",NumEvents," events using ",NumPseudoExp," pseudo-experiments"
!      write(*,"(A,I2,A,I6,A,I9,A)") "reading histogram ",iHisto," for likelihood analysis with Lumi=",Lumi,"fb^-1 using ",NumPseudoExp," pseudo-experiments"
      write(*,"(A,I2,A,I6,A,I9,A,F9.6)") "reading histogram ",iHisto," for likelihood analysis with Lumi=",Lumi,"fb^-1 using ",NumPseudoExp," pseudo-experiments, and scale uncertainty", DeltaN
      write(14,"(A,I2,A,I6,A,I9,A,F9.6)") "# reading histogram ",iHisto," for likelihood analysis with Lumi=",Lumi,"fb^-1 using ",NumPseudoExp," pseudo-experiments, and scale uncertainty", DeltaN
      write(*,*) "Writing output to file : ", filename_out
      write(*,*) ""
!       call likelihood(iHisto,NumEvents,NumPseudoExp)
      call likelihood(iHisto,Lumi,NumPseudoExp,DeltaN)

  elseif( op.eq.9 ) then! SYNTAX: NHi dV1 dA1 file1 ... dV6 dA6 file6 outfilenametag
      VcouplSM = 0.24364d0
      AcouplSM = 0.60069d0
      write(*,"(2X,A,F8.5,A,F8.5)") "Fixed SM t-tb-Z couplings V=",VcouplSM," A=",AcouplSM
      call GetArg(2,iHisto_str)
      read(iHisto_str,"(I2)") iHisto   
      do NArg=3,NumArgs-1,3
            call GetArg(NArg,dVcoupl_str)
            read(dVcoupl_str,"(F5.2)") dVcoupl(NArg/3)
            call GetArg(NArg+1,dAcoupl_str)
            read(dAcoupl_str,"(F5.2)") dAcoupl(NArg/3)
            call GetArg(NArg+2,filename_in)
            open(unit=10+(NArg/3),file=trim(filename_in),form='formatted',access='sequential')! units=11,...,16
            write(*,"(2X,A,I1,X,A,A,F6.3,A,F6.3)") "Reading file ",(NArg/3),trim(filename_in), " with dV=",dVcoupl(NArg/3)," dA=",dAcoupl(NArg/3)
      enddo
      call GetArg(NumArgs,filename_out)
      open(unit=10,file=trim(filename_out),form='formatted',access='sequential',status='replace') ! open output file (last one in arguments list)
      call ttbzcoupl(iHisto,VcouplSM,AcouplSM,dVcoupl,dAcoupl,filename_out)
      close(11);close(12);close(13);close(14);close(15);close(16)
  endif


! close all files
  close(10)


CONTAINS 




SUBROUTINE add_avg(NumArgs,op,silent)
implicit none
integer :: NArg,row,NumArgs,op
character(len=*),parameter :: fmt1 = "(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)"
integer,parameter :: MaxFiles=50
integer :: NHisto(1:MaxFiles)=-999999,Hits(1:MaxFiles)=-999999
real(8) :: BinVal(1:MaxFiles)=-1d-99,Value(1:MaxFiles)=-1d-99,Error(1:MaxFiles)=-1d-99
real(8) :: SumValue,SumError,n
integer :: SumHits
character :: dummy*(1)
logical :: silent

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
    if( .not.silent ) write(* ,fmt1) NHisto(1),"|",BinVal(1),"|",SumValue,"|",SumError,"|",SumHits,"|"
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
         write(*  ,fmt1) NHisto(1),"|",BinVal(1),"|",newvalue,"|",newerror,"|",newhits,"|"

!          do ibin = 1, nrebin
!             write(13 ,fmt1) NHisto(1),"|",BinVal(1+ibin),"|",newvalue,"|",newerror,"|",newhits,"|"
!          enddo

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





SUBROUTINE likelihood(iHisto,LumiORNumEvents,NumPseudoExp,DeltaN)
implicit none
integer :: iHisto,LumiORNumEvents,NumPseudoExp
real(8) :: DeltaN
integer,parameter :: MaxBins=1000, MaxEvents=200000
integer,parameter :: NumLLbins=1000! the LLratio should only vary between - and + this number
character(len=*),parameter :: fmt1 = "(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)"
character :: dummy*(1)
integer :: iBin,iPseudoExp,NumBins1=0,NumBins2=0,ranBin,BinAccept(1:MaxEvents),NumAccepted,NEvent,TryEvents
integer :: NumEvents,NumExpectedEvents(1:2),iHypothesis
integer :: NHisto(1:2)=-999999,Hits(1:2,1:MaxBins)=-999999,WhichBin
real(8) :: Lumi,BinVal(1:2,1:MaxBins)=-1d-99,Value(1:2,1:MaxBins)=-1d-99,Error(1:2,1:MaxBins)=-1d-99
real(8) :: ymax,zran(1:2),nran(1:2),ValAccepted(1:MaxEvents),LLRatio,BinSize,MaxPoisson(1:2)
real(8) :: sigmatot(1:2),check(1:2),alpha(1:2,1:MaxBins),IntLLRatio(1:2),checkalpha(1:MaxBins)
integer :: SelectedEvent,LLbin,PlotObsEvts(1:2,1:10000),i,s,offset,SUA
real(8) :: alphamin,alphamax,betamin,betamax,rescale(1:2),sran
real(8) :: LLRatio_array(1:2,1:NumPseudoExp),LLRatio_min,LLRatio_max
logical ::  GotNumEvents,PoissonEvents,PoissonBins,useshape
type :: Histogram
    integer :: NBins
    real(8) :: BinSize
    real(8) :: LowVal
    integer :: Hits(1:1000)
end type
type(Histogram) :: LLHisto(1:2)
logical, parameter :: normalize=.true.



! NumEvents=LumiORNumEvents
Lumi = dble(LumiORNumEvents)


!--------------------------------------------------
!          0. init
!--------------------------------------------------
LLHisto(1)%NBins   = 1000
LLHisto(1)%BinSize = 0.5d0
LLHisto(1)%LowVal  = -200d0
LLHisto(1)%Hits(:) = 0

LLHisto(2)%NBins   = LLHisto(1)%NBins
LLHisto(2)%BinSize = LLHisto(1)%BinSize
LLHisto(2)%LowVal  = LLHisto(1)%LowVal
LLHisto(2)%Hits(:) = LLHisto(1)%Hits(:)
PlotObsEvts=0

PoissonEvents=.true.
!  Scale Uncertainty Approach: 
! 1=rescale PREDICTED cross-secs and distr (conservative)
! 2=change OBSERVED events in each pseudoexp
SUA=1                   

if (SUA .ne. 1 .and. SUA .ne. 2 .and. DeltaN .ne. 1d0) then
   print *, "WARNING : SCALE UNCERTAINTY NOT USED!"
   stop
endif



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
  write(*,"(A,1PE16.8,A,I6)") "Total cross section of input file 1: ",sigmatot(1),"   <-->   Number of events: ",NumExpectedEvents(1)
  write(*,"(A,1PE16.8,A,I6)") "Total cross section of input file 2: ",sigmatot(2),"   <-->   Number of events: ",NumExpectedEvents(2)
  check(1:2) = 0d0

  rescale=1d0
  if (SUA .eq. 1) then
     ! Added -- scale uncertainty -- this should also change 
     if (sigmatot(1) .gt. sigmatot(2)) then
        if (sigmatot(1)/DeltaN .gt. sigmatot(2)*DeltaN) then
           sigmatot(1)=sigmatot(1)/DeltaN
           sigmatot(2)=sigmatot(2)*DeltaN
           rescale(1)=1d0/DeltaN
           rescale(2)=DeltaN
        else
           rescale(1)=1d0
           rescale(2)=sigmatot(1)/sigmatot(2)
           sigmatot(2)=sigmatot(1)
        endif
     else
        if (sigmatot(1)*DeltaN .lt. sigmatot(2)/DeltaN) then
           sigmatot(1)=sigmatot(1)*DeltaN
           sigmatot(2)=sigmatot(2)/DeltaN
           rescale(1)=DeltaN
           rescale(2)=1d0/DeltaN
        else
           rescale(1)=1d0
           rescale(2)=sigmatot(1)/sigmatot(2)
           sigmatot(2)=sigmatot(1)
        endif
     endif
  endif

  NumExpectedEvents(1) = int(Lumi * sigmatot(1) *8 ) ! factor of 8 from lepton species
  NumExpectedEvents(2) = int(Lumi * sigmatot(2) *8 ) 


  write(*,"(A,I6,A)") "Scale uncertainty for null hypothesis-> ", NumExpectedEvents(1) , " events"
  write(*,"(A,I6,A)") "Scale uncertainty for alt hypothesis-> ", NumExpectedEvents(2) , " events"
  
  MaxPoisson(1) = Poisson(NumExpectedEvents(1),NumExpectedEvents(1))
  MaxPoisson(2) = Poisson(NumExpectedEvents(2),NumExpectedEvents(2))

  print *, 'Max of Poisson:', MaxPoisson(1:2)
! normalize distributions
  if( normalize ) then
   print *, "Normalizing input files" 
   do iBin=1,NumBins1! normalize histogram
       Value(1,iBin) = Value(1,iBin)*BinSize/sigmatot(1)*rescale(1)
       Value(2,iBin) = Value(2,iBin)*BinSize/sigmatot(2)*rescale(2)
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
  call random_seed()     
do iHypothesis=1,2

   if (.not. PoissonEvents) then       ! dont incl cross-section info, so same number of events for every pseudoexperiment
      !   NumEvents = NumExpectedEvents(1)              ! choose number of events fixed for each hypothesis
      NumEvents = NumExpectedEvents(iHypothesis)    ! choose number of events differently for each hypothesis

      write(*,"(A,I6,A)") "Generating ",NumEvents," events"
      if( NumEvents.gt.MaxEvents ) then
         print *, "Error: NumEvents is too large. Increase MaxEvents in SUBROUTINE likelihood."
         stop
      endif
   else
      write(*,*)  "Generating Poisson distributed events for each pseudoexperiment"
   endif

   do iPseudoExp=1,NumPseudoExp

      if( mod(iPseudoExp,10000).eq.0 ) print *, "Pseudo experiment ",iPseudoExp,"/",NumPseudoExp

      
      if (PoissonEvents) then
               
!*************************************************************
!  2.0 generate number of events for null hypothesis
!*************************************************************
         GotNumEvents=.false.
               
         do while (.not. GotNumEvents) 
            call random_number(nran(1:2))
            nran(1)=nran(1)*MaxPoisson(iHypothesis)
            TryEvents=int(2d0*NumExpectedEvents(iHypothesis)*nran(2)) 
            if ( Poisson(NumExpectedEvents(iHypothesis),TryEvents) .gt. nran(1) ) then
               NumEvents = TryEvents
               GotNumEvents=.true.
!               PlotObsEvts(iHypothesis,TryEvents)=PlotObsEvts(iHypothesis,TryEvents)+1     ! this is to check the Poisson dist.
            endif
         enddo
         if( NumEvents.gt.MaxEvents ) then
            print *, "Error: NumEvents is too large. Increase MaxEvents in SUBROUTINE likelihood."
            stop
         endif
      endif

      if (SUA .eq. 2) then  
         !************************************************************
         !  2.0a change number of observed events for scale uncertainty
         !************************************************************    
         
         call random_number(sran)
         sran=sran*(DeltaN**2-1d0)/DeltaN + 1d0/DeltaN
         write(201,*) sran
         ! this should be uniformly distr in [1/DeltaN,DeltaN] -- check
         NumEvents=NumEvents*sran
         PlotObsEvts(iHypothesis,NumEvents)=PlotObsEvts(iHypothesis,NumEvents)+1     ! this is to check the Poisson dist.
      endif



!************************************************************
!  2.1 generate event sample according to hypothesis
!************************************************************
      ! find max. of y-axis
      ymax=-1d13;
      do iBin=1,NumBins1
         if(  Value(iHypothesis,iBin).gt.ymax ) ymax = Value(iHypothesis,iBin)
      enddo

      NumAccepted=0
!      call random_seed()     
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
!      LLRatio = 0d0
       if (PoissonEvents) then
         LLRatio=NumEvents*dlog(1d0*NumExpectedEvents(1)/NumExpectedEvents(2))
      else
         LLRatio = 0d0
      endif
      do NEvent=1,NumEvents
         SelectedEvent = BinAccept(NEvent)
         if( ValAccepted(NEvent).ne.Value(iHypothesis,SelectedEvent) ) then! can be removed/simplified later 
             print *, "error in ll"
             stop
         endif
         LLRatio = LLRatio + 2d0*dlog( (Value(1,SelectedEvent))/(Value(2,SelectedEvent)) )
      enddo
      if (iHypothesis .eq. 1 .and. iPseudoExp .eq. 1) then
         offset=-100d0*(int(LLRatio)/100)

      endif
      LLRatio=LLRatio+offset
      if (iPseudoExp .eq. 1) then
         print *, 'LLRatio=',LLRatio,offset
      endif


      LLRatio_array(iHypothesis,iPseudoExp)=LLRatio
  enddo! iPseudoExp
  do i=1,2000
     s=101+iHypothesis
     write(s,*) i,PlotObsEvts(iHypothesis,i)
  enddo

enddo! iHypothesis

LLRatio_max=-1d-6
LLRatio_min=1d6

do iHypothesis=1,2
   do iPseudoExp=1,NumPseudoExp
      if (LLRatio_array(iHypothesis,iPseudoExp) .gt. LLRatio_max) then
         LLRatio_max=LLRatio_array(iHypothesis,iPseudoExp) 
      endif
      if (LLRatio_array(iHypothesis,iPseudoExp) .lt. LLRatio_min) then
         LLRatio_min=LLRatio_array(iHypothesis,iPseudoExp) 
      endif
   enddo
enddo
print *, LLRatio_min,LLRatio_max
!LLRatio_max=100d0*(int(LLRatio_max/100)+1)
!LLRatio_min=100d0*(int(LLRatio_max/100)-1)
LLRatio_max=ceiling(LLRatio_max)
LLRatio_min=floor(LLRatio_min)

if (LLRatio_max .lt. LLRatio_min) then
   print *, 'ERROR: MAX OF LL RATIO IS SMALLER THAN MIN!'
   print *, LLRatio_max,LLRatio_min
   stop
endif


LLHisto(1)%NBins   = 1000
LLHisto(1)%BinSize = (LLRatio_max-LLRatio_min)/LLHisto(1)%NBins
LLHisto(1)%LowVal  = LLRatio_min
LLHisto(1)%Hits(:) = 0
print *, LLHisto(1)%BinSize
print *, LLRatio_min, LLRatio_max
pause
LLHisto(2)%NBins   = LLHisto(1)%NBins
LLHisto(2)%BinSize = LLHisto(1)%BinSize
LLHisto(2)%LowVal  = LLHisto(1)%LowVal
LLHisto(2)%Hits(:) = LLHisto(1)%Hits(:)

!************************************************************
!  2.3 bin the likelihood value 
!************************************************************
do iHypothesis=1,2
   do iPseudoExp=1,NumPseudoExp
      LLRatio=LLRatio_array(iHypothesis,iPseudoExp)
      WhichBin = (LLRatio-LLHisto(iHypothesis)%LowVal)/LLHisto(iHypothesis)%BinSize + 1
      WhichBin=int(WhichBin)
      if( WhichBin.lt.0 ) WhichBin = 1
      if( WhichBin.gt.LLHisto(iHypothesis)%NBins ) WhichBin = LLHisto(iHypothesis)%NBins
      LLHisto(iHypothesis)%Hits(WhichBin) = LLHisto(iHypothesis)%Hits(WhichBin) + 1
   enddo

         
!************************************************************
! 2.4 Find the integral under the LL distribution
!************************************************************

      do LLbin=1,LLHisto(iHypothesis)%NBins
         IntLLRatio(iHypothesis)=IntLLRatio(iHypothesis)+&
              LLHisto(iHypothesis)%BinSize * LLHisto(iHypothesis)%Hits(LLbin)
         alpha(iHypothesis,LLBin)=IntLLRatio(iHypothesis)/(NumPseudoExp*LLHisto(iHypothesis)%BinSize)
      enddo
enddo
!************************************************************
!  3. write out the LL distribution, and find integral under each distribution
!************************************************************
      print *, ""
      print *, "Writing log-likelihood distribution to outful file"
      print *, ""
      do LLbin=1,LLHisto(1)%NBins
!         print *, LLBin
           write(14,"(2X,I4,2X,1PE16.8,2X,I10,2X,1PE16.8,2X,I10,2X,1PE16.8,2X,1PE16.8)") LLbin, LLHisto(1)%LowVal+LLbin*LLHisto(1)%BinSize, LLHisto(1)%Hits(LLbin),LLHisto(2)%LowVal+LLbin*LLHisto(2)%BinSize, LLHisto(2)%Hits(LLbin),alpha(1,LLBin),alpha(2,LLBin)
      enddo


!************************************************************
!  4. Now find the point at which alpha=1-beta
!************************************************************

! for the time being, I'm going to assume that both distributions have the same range and binning - can change this later
   do LLbin=1,LLHisto(1)%NBins
      checkalpha(LLBin)=alpha(1,LLBin)+alpha(2,LLBin)-1d0
      print *, LLBIn,alpha(1,LLBin),alpha(2,LLBin),checkalpha(LLBin)
      if (checkalpha(LLBin)*checkalpha(LLBin-1) .lt. 0d0) then
         alphamin=alpha(1,LLBin-1)
         alphamax=alpha(1,LLBin)
         betamin =alpha(2,LLBin-1)
         betamax =alpha(2,LLBin)
         
         if (checkalpha(LLBin) .eq. 1d0 .and. checkalpha(LLBin-1) .eq.-1d0) then
            ! this is comparing two identical hypotheses!
            alphamin=1d0
            alphamax=1d0
            betamin=0d0
            betamax=0d0
         endif
      endif
   enddo

   print *, 'alpha value in range:', alphamin,alphamax
   print *, 'betaa value in range:', betamin,betamax 

   write(14,"(A,2X,1PE16.8,2X,1PE16.8)") "# alpha value in range:",alphamin,alphamax
   write(14,"(A,2X,1PE16.8,2X,1PE16.8)") "# beta value in range:",betamin,betamax


END SUBROUTINE






SUBROUTINE ttbzcoupl(iHisto,VcouplSM,AcouplSM,dVcoupl,dAcoupl,filename_out)
IMPLICIT NONE
integer :: iHisto
real(8) :: VcouplSM,AcouplSM,dVcoupl(1:6),dAcoupl(1:6),Vcoupl,Acoupl
character :: filename_out*(50),NewFileName*(100)
character(len=*),parameter :: fmt1 = "(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)"
character :: dummy*(1)
integer,parameter :: MaxBins=100
integer :: NHisto=-999999,Hits=-999999,TheUnit,ibin,NumBins,Vgrid,Agrid,Ai,Vi
real(8) :: BinVal(1:6,0:MaxBins)=-1d-99,Value(1:6,0:MaxBins)=-1d-99,Error=-1d-99
real(8) :: res,CoefFits(1:6),couplMatrix_new(1:6),dV_new,dA_new,NewValue,dVrange,dArange,dVstep,dAstep
complex(8) :: couplMatrix(1:6,1:7)

  ! reading all input files
  do TheUnit=11,16! looping over all 6 input files
      ibin=0
      do while(.not.eof(TheUnit))! looping over all lines in the input file
          read(unit=TheUnit,fmt="(A)") dummy; if(dummy(1:1).eq."#") cycle; backspace(unit=TheUnit) ! skip # lines
          read(unit=TheUnit,fmt="(I2)") NHisto
          if(NHisto.eq.iHisto) then
              backspace(unit=TheUnit)
              ibin=ibin+1
              read(unit=TheUnit,fmt=fmt1) NHisto,dummy,BinVal(TheUnit-10,ibin),dummy,Value(TheUnit-10,ibin),dummy,Error,dummy,Hits,dummy
          endif
      enddo 
  enddo
  NumBins=ibin



dummy=" "
! write(*,"(2X,A,") ""
do ibin=1,NumBins! looping over all histogram bins

! Solving the system of 6 linear equations for coefficients of A^2, V^2, 1, A*V, A, V
  do TheUnit=1,6
     if( BinVal(TheUnit,ibin).ne.BinVal( mod(TheUnit,6)+1,ibin)  ) print *, "ERROR: Different bin values!" 
     Vcoupl = VcouplSM * ( 1d0 + dVcoupl(TheUnit) )
     Acoupl = AcouplSM * ( 1d0 + dAcoupl(TheUnit) )
     res = Value(TheUnit,ibin)
     couplMatrix(TheUnit,1:7) = (/ dcmplx(Acoupl**2), dcmplx(Vcoupl**2), dcmplx(1d0), dcmplx(Acoupl*Vcoupl), dcmplx(Acoupl), dcmplx(Vcoupl), dcmplx(res) /)
!      write(*,"(2X,A,6F8.5,F10.6)") "Matrix ",dble(couplMatrix(TheUnit,1:7))
   enddo
   CoefFits(1:6) = dble( go_Gauss_64(6,couplMatrix(1:6,1:7)) )
!  CoefFits(1:6) = dble( go_GaussLU(6,couplMatrix(1:6,1:7)) )
   write(*,"(2X,A,I2,A,6F10.6)") "Bin=",ibin,". Fit coeffs.:",CoefFits(1:6)


!--- GRID SETUP -----------
     dVrange=3.0d0   ! this means variation by +/- dVrange
     dArange=0.4d0   ! this means variation by +/- dArange
     Vgrid=60
     Agrid=40
!--------------------------

     dVstep=dVrange/dble(Vgrid)
     dAstep=dArange/dble(Agrid)
     dV_new=0d0; dA_new=0d0
     do Vi=0,Vgrid
        dV_new = Vi * dVstep
        do Ai=0,Agrid
           dA_new = Ai * dAstep


!           generating the +V +A results        
            Vcoupl = VcouplSM * ( 1d0 + dV_new )
            Acoupl = AcouplSM * ( 1d0 + dA_new )
            couplMatrix_new(1:6) = (/ dcmplx(Acoupl**2), dcmplx(Vcoupl**2), dcmplx(1d0), dcmplx(Acoupl*Vcoupl), dcmplx(Acoupl), dcmplx(Vcoupl) /)
            NewValue = couplMatrix_new(1)*CoefFits(1) + couplMatrix_new(2)*CoefFits(2) + couplMatrix_new(3)*CoefFits(3)  &
                     + couplMatrix_new(4)*CoefFits(4) + couplMatrix_new(5)*CoefFits(5) + couplMatrix_new(6)*CoefFits(6)
            write(NewFileName,"(A,A,I2,A,1F4.2,A,1F4.2,A)") trim(filename_out),"_Hi",iHisto,"_V+",dV_new,"_A+",dA_new,".dat"
            if( ibin.eq.1 ) then
                open(unit=10,file=trim(NewFileName),form='formatted',access='sequential')
            else
                open(unit=10,file=trim(NewFileName),form='formatted',access='sequential',position='append')
            endif
            write(unit=10,fmt=fmt1) iHisto,dummy,BinVal(1,ibin),dummy,NewValue,dummy,0d0,dummy,0,dummy
            close(10)


!           generating the +V -A results        
            Vcoupl = VcouplSM * ( 1d0 + dV_new )
            Acoupl = AcouplSM * ( 1d0 - dA_new )
            couplMatrix_new(1:6) = (/ dcmplx(Acoupl**2), dcmplx(Vcoupl**2), dcmplx(1d0), dcmplx(Acoupl*Vcoupl), dcmplx(Acoupl), dcmplx(Vcoupl) /)
            NewValue = couplMatrix_new(1)*CoefFits(1) + couplMatrix_new(2)*CoefFits(2) + couplMatrix_new(3)*CoefFits(3)  &
                     + couplMatrix_new(4)*CoefFits(4) + couplMatrix_new(5)*CoefFits(5) + couplMatrix_new(6)*CoefFits(6)
            write(NewFileName,"(A,A,I2,A,1F4.2,A,1F4.2,A)") trim(filename_out),"_Hi",iHisto,"_V+",dV_new,"_A-",dA_new,".dat"
            if( ibin.eq.1 ) then
                open(unit=10,file=trim(NewFileName),form='formatted',access='sequential')
            else
                open(unit=10,file=trim(NewFileName),form='formatted',access='sequential',position='append')
            endif
            write(unit=10,fmt=fmt1) iHisto,dummy,BinVal(1,ibin),dummy,NewValue,dummy,0d0,dummy,0,dummy
            close(10)


!           generating the -V +A results        
            Vcoupl = VcouplSM * ( 1d0 - dV_new )
            Acoupl = AcouplSM * ( 1d0 + dA_new )
            couplMatrix_new(1:6) = (/ dcmplx(Acoupl**2), dcmplx(Vcoupl**2), dcmplx(1d0), dcmplx(Acoupl*Vcoupl), dcmplx(Acoupl), dcmplx(Vcoupl) /)
            NewValue = couplMatrix_new(1)*CoefFits(1) + couplMatrix_new(2)*CoefFits(2) + couplMatrix_new(3)*CoefFits(3)  &
                     + couplMatrix_new(4)*CoefFits(4) + couplMatrix_new(5)*CoefFits(5) + couplMatrix_new(6)*CoefFits(6)
            write(NewFileName,"(A,A,I2,A,1F4.2,A,1F4.2,A)") trim(filename_out),"_Hi",iHisto,"_V-",dV_new,"_A+",dA_new,".dat"
            if( ibin.eq.1 ) then
                open(unit=10,file=trim(NewFileName),form='formatted',access='sequential')
            else
                open(unit=10,file=trim(NewFileName),form='formatted',access='sequential',position='append')
            endif
            write(unit=10,fmt=fmt1) iHisto,dummy,BinVal(1,ibin),dummy,NewValue,dummy,0d0,dummy,0,dummy
            close(10)


!           generating the -V -A results        
            Vcoupl = VcouplSM * ( 1d0 - dV_new )
            Acoupl = AcouplSM * ( 1d0 - dA_new )
            couplMatrix_new(1:6) = (/ dcmplx(Acoupl**2), dcmplx(Vcoupl**2), dcmplx(1d0), dcmplx(Acoupl*Vcoupl), dcmplx(Acoupl), dcmplx(Vcoupl) /)
            NewValue = couplMatrix_new(1)*CoefFits(1) + couplMatrix_new(2)*CoefFits(2) + couplMatrix_new(3)*CoefFits(3)  &
                     + couplMatrix_new(4)*CoefFits(4) + couplMatrix_new(5)*CoefFits(5) + couplMatrix_new(6)*CoefFits(6)
            write(NewFileName,"(A,A,I2,A,1F4.2,A,1F4.2,A)") trim(filename_out),"_Hi",iHisto,"_V-",dV_new,"_A-",dA_new,".dat"
            if( ibin.eq.1 ) then
                open(unit=10,file=trim(NewFileName),form='formatted',access='sequential')
            else
                open(unit=10,file=trim(NewFileName),form='formatted',access='sequential',position='append')
            endif
            write(unit=10,fmt=fmt1) iHisto,dummy,BinVal(1,ibin),dummy,NewValue,dummy,0d0,dummy,0,dummy
            close(10)
        enddo
     enddo

enddo! ibin


END SUBROUTINE




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










FUNCTION go_Gauss_64(N,matrix)
implicit none
complex(8), intent(in) :: matrix(1:N,1:N+1)
integer, intent(in) :: N
complex(8) :: go_Gauss_64(1:N)
complex(8) :: work(1:N,1:N+1)
integer :: i,j
complex(8) :: mult, tsum

    work(1:N,1:N+1) = matrix(1:N,1:N+1)
    do i = 1,N
      if (cdabs(work(i,i)) <= 1.0e-6) then
         print *, i,N
         print *, "go_Gauss, Zero pivot element: ",cdabs(work(i,i))
      endif
      do j = i+1,N
        mult = work(j,i)/work(i,i)
        work(j,1:N+1) = work(j,1:N+1) - mult*work(i,1:N+1)
      end do
    end do
    do i=N,1,-1
      tsum = work(i,N+1)
      do j=i+1,N
        tsum = tsum - work(i,j)*work(j,N+1)
      end do
      work(i,N+1) = tsum/work(i,i)
    end do
    go_Gauss_64(1:N) = work(1:N,N+1)
    return
END FUNCTION








END PROGRAM












