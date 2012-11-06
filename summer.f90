PROGRAM summer
use ifport
implicit none
integer :: NumArgs,NArg,op,iHisto
character :: filename_in*(50),filename_out*(50),operation*(10),factor_str*(20),iHisto_str*(5)
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
  else
      print *, "operation not available (add,avg,sumbins,mul,quo,rebin)",trim(operation)
      stop
  endif

! check for invalid arguments
  if( (op.lt.3 .and. NumArgs.le.3) .or. (op.eq.3 .and. NumArgs.ne.3) .or. (op.eq.4 .and. NumArgs.ne.4) .or. (op.eq.5 .and. NumArgs.le.3) .or. NumArgs.gt.49) then
    print *, "invalid number of arguments",op,NumArgs
    print *, "the first argument specifies the operation, the last argument is the output file"
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
      call rebinning(iHisto)
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



SUBROUTINE rebinning(iHisto)
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


      read(unit=12,fmt="(A)") dummy
      if(dummy(1:1).eq."#") cycle
      backspace(unit=12) ! go to the beginning of the line
      read(unit=12,fmt=fmt1) NHisto(1),dummy,BinVal(1),dummy,Value(1),dummy,Error(1),dummy,Hits(1),dummy

      if(NHisto(1).eq.iHisto) then
!         write(* ,fmt1) NHisto(1),"|",BinVal(1),"|",Value(1),"|",Error(1),"|",Hits(1),"|"

         read(unit=12,fmt="(A)") dummy
         if(dummy(1:1).eq."#") cycle
         backspace(unit=12) ! go to the beginning of the line

         read(unit=12,fmt=fmt1) NHisto(2),dummy,BinVal(2),dummy,Value(2),dummy,Error(2),dummy,Hits(2),dummy
!         write(* ,fmt1) NHisto(2),"|",BinVal(2),"|",Value(2),"|",Error(2),"|",Hits(2),"|"
    
         write(13 ,fmt1) NHisto(1),"|",BinVal(1),"|",(Value(1)+Value(2))/2d0,"|",dsqrt(Error(1)**2+Error(2)**2)/2d0,"|",Hits(1)+Hits(2),"|"
         write(*  ,fmt1) NHisto(1),"|",BinVal(1),"|",(Value(1)+Value(2))/2d0,"|",dsqrt(Error(1)**2+Error(2)**2)/2d0,"|",Hits(1)+Hits(2),"|"
      else
         write(13 ,fmt1) NHisto(1),"|",BinVal(1),"|",Value(1),"|",Error(1),"|",Hits(1),"|"
         write(*  ,fmt1) NHisto(1),"|",BinVal(1),"|",Value(1),"|",Error(1),"|",Hits(1),"|"
      endif

  enddo


END SUBROUTINE



