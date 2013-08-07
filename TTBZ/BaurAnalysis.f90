program BaurAnalysis
! RR 06 Aug 2013
! to get P-value, the following arguments are used (in order) : 
! 'GetPVal' -- tells the program what to do
! SM_infile  (string) -- input file with SM distributions
! BSM_infile (string) -- input file with 'BSM' (i.e. using non-standard values for F_1,2_V,A) distributions
! Histo (int) -- Label of observable in SM_infile and BSM_infile to be used to obtain P-value
! Integrated luminosity (real)
! DeltaN (real) -- scale (or other) uncertainty  
  use ifport
  implicit none
  integer,parameter :: MaxBins=200
  integer :: NumArgs,Histo,i,j,k,NumBins_SM,NumBins_BSM,NumBins,sum_Events_SM,NumEventBins
  integer :: NHisto(1:MaxBins),Hits(1:MaxBins),Events_SM(1:MaxBins),Events_BSM(1:MaxBins)
  real(8) :: PreFac,data,DeltaN,BinSize_SM,BinSize_BSM,BinSize,fbar_sq,f,chisq,P
  real(8) :: BinVal(1:MaxBins),Valu_SM(1:MaxBins),Error_SM(1:MaxBins),Valu_BSM(1:MaxBins),Error_BSM(1:MaxBins)
  character ::  operation*(10), SM_infile*(50), BSM_infile*(50),dummy*(1),Histo_str*(5),data_str*(5),DeltaN_str*(5)
character(len=*),parameter :: fmt1 = "(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)"



! get number of arguments
  NumArgs = NArgs()-1

! Prefactor to multiply all the data
  PreFac = 2d0*2d0*2d0   ! 2 lepton families from Z decay, two lepton families from t decay, 2 to allow other top to decay hadronically
  
! to recover the Baur analysis  - tagging efficiencies
  PreFac = PreFac/2d0*0.85d0**3*0.4d0
  
  call GetArg(1,operation)

  if( trim(operation) .ne. 'GetPVal') then
     print *, 'Operation has to be GetPVal'
     stop
  endif

  if (NumArgs .ne. 6)   then
     print *, 'Requires six arguments, you entered ', NumArgs
     stop
  endif

  call GetArg(2,SM_infile)
  call GetArg(3,BSM_infile)
  call GetArg(4,Histo_str)
  call GetArg(5,data_str)
  call GetArg(6,DeltaN_str)
  read(Histo_str,"(I2)") Histo
  read(data_str,"(F10.6)") data
  read(DeltaN_str,"(F10.6)") DeltaN
  
  PreFac=PreFac*data

  open(unit=11,file=trim(SM_infile),form='formatted',access='sequential')
  open(unit=12,file=trim(BSM_infile),form='formatted',access='sequential')

  i=0
  do while(.not.eof(11))  ! loop over all rows

      read(unit=11,fmt="(A)") dummy
      if(dummy(1:1).eq."#") cycle
      backspace(unit=11) ! go to the beginning of the line
      i=i+1
      read(unit=11,fmt=fmt1) NHisto(i),dummy,BinVal(i),dummy,Valu_SM(i),dummy,Error_SM(i),dummy,Hits(i),dummy

      if(NHisto(i).ne.Histo) then
         i=i-1
      endif
   enddo
   NumBins_SM=i
   
   BinSize_SM=BinVal(2)-BinVal(1)

   i=0
   do while(.not.eof(12))  ! loop over all rows

      read(unit=12,fmt="(A)") dummy
      if(dummy(1:1).eq."#") cycle
      backspace(unit=12) ! go to the beginning of the line
      i=i+1
      read(unit=12,fmt=fmt1) NHisto(i),dummy,BinVal(i),dummy,Valu_BSM(i),dummy,Error_BSM(i),dummy,Hits(i),dummy

      if(NHisto(i).ne.Histo) then
         i=i-1
      endif
   enddo
   NumBins_BSM=i

   BinSize_BSM=BinVal(2)-BinVal(1)
   if (BinSize_BSM .ne. BinSize_SM) then
      print *, 'Error: bins sizes must be same for the two input files!'
      stop
   endif
   if (NumBins_SM .ne. NumBins_BSM) then
      print *, 'Error: number of bins must be the same for the two input files!'
      print *, NumBins_SM,  NumBins_BSM
      stop
   endif

   NumBins=NumBins_SM
   Valu_SM  = Valu_SM  * PreFac * BinSize_SM
   Valu_BSM = Valu_BSM * PreFac * BinSize_SM
   
   j=0
   sum_Events_SM=0
   fbar_sq=0d0
   
   do i = 1,NumBins
      if (Valu_SM(i) .ge. 5d0 .and. Valu_BSM(i) .ge. 5d0) then
!      if (Valu_SM(i) .ge. 1d0 .and. Valu_BSM(i) .ge. 1d0) then
         j=j+1
         Events_SM(j) =floor(Valu_SM(i))
         Events_BSM(j)=floor(Valu_BSM(i))
         sum_Events_SM=sum_Events_SM+Events_SM(j)
         fbar_sq = fbar_sq + 1.d0*Events_BSM(j)**2/Events_SM(j)
      endif
   enddo
   NumEventBins=j
   fbar_sq=fbar_sq/sum_Events_SM


   if ( dsqrt(fbar_sq) .gt. 1d0+DeltaN ) then
      f=1d0+DeltaN
   elseif ( dsqrt(fbar_sq) .lt. 1d0/(1d0+DeltaN) ) then
      f=1d0/(1d0+DeltaN)
   else
      f=dsqrt(fbar_sq)
   endif
   
   print *, 'fbar_sq', fbar_sq
   print *, 'f = ', f
!! RR unsure about this
   chisq=NumEventBins-1
!!   chisq=0
!!!
   do j=1,NumEventBins
      chisq=chisq + (Events_BSM(j)*1.0d0 - f*Events_SM(j))**2/(f*Events_SM(j))
      print *, Events_BSM(j),Events_SM(j)
   enddo
   
!! RR UNSURE OF THIS
!!   chisq=chisq/(NumEventBins-1)
!!!
!!   chisq=50d0
   print *, 'chisq = ', chisq
   if (mod(NumEventBins,2).ne.0) then
      print *, 'WARNING: ODD NUMBER OF BINS'
      print *, 'The expression for the P-value is no longer correct'
   endif
   ! now that I have chi-sq, I can get p-value
   P=0d0
   do k=0,NumEventBins/2-1
      P=P+(chisq/2d0)**k/factorial(k)*exp(-chisq/2d0)
   enddo
   print *, 'P-value: ', P
   



contains 
  FUNCTION factorial(N)
    implicit none
!    integer Factorial,n,i
    integer n,i
    real(8) Factorial
    
    factorial = 1d0
    do i=n,1,-1
       factorial = factorial * i
    end do
    return
  END FUNCTION factorial

 end program BaurAnalysis

      
   
      


