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
  integer :: NumArgs,Histo1,Histo2,NumBins_SM,NumBins_BSM,NumBins,NumEventBins
  integer :: NumBins_SM2,NumBins_BSM2,NumBins2,NumEventBins2,i
  integer :: Events_SM(1:MaxBins),Events_BSM(1:MaxBins)
  integer :: Events_SM2(1:MaxBins),Events_BSM2(1:MaxBins)
  real(8) :: PreFac,data,DeltaN,BinSize_SM,chisq,P,sigs
  real(8) :: Valu_SM(1:MaxBins),Valu_BSM(1:MaxBins)
  real(8) :: Valu_SM2(1:MaxBins),Valu_BSM2(1:MaxBins)
  character :: DelF1V*(10), DelF1A*(10),cachestr*(10)
  integer :: cache
  character ::  operation*(10), SM_infile*(50), BSM_infile*(50),dummy*(1),Histo1_str*(5),Histo2_str*(5),data_str*(5),DeltaN_str*(5)
character(len=*),parameter :: fmt1 = "(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)"



! get number of arguments
  NumArgs = NArgs()-1

! Prefactor to multiply all the data
  PreFac = 2d0*2d0*2d0   ! 2 lepton families from Z decay, two lepton families from t decay, 2 to allow other top to decay hadronically
  
! to recover the Baur analysis  - tagging efficiencies
  PreFac = PreFac/2d0*0.85d0**3*0.4d0
  
  call GetArg(1,operation)

  if( (trim(operation) .ne. 'GetPVal') .and. (trim(operation) .ne. 'GetLogL') ) then
     print *, 'Operation has to be GetPVal or GetLogL'
     stop
  endif

  if (NumArgs .lt. 6)   then
     print *, 'Requires six arguments, you entered ', NumArgs
     stop
  endif

  call GetArg(2,SM_infile)
  call GetArg(3,BSM_infile)
  call GetArg(4,Histo1_str)
  call GetArg(5,Histo2_str)
  call GetArg(6,data_str)
  call GetArg(7,DeltaN_str)
  call GetArg(8,Delf1V)     ! these are just for output, not needed??
  call GetArg(9,Delf1A)
  call GetArg(10,cachestr)
  read(Histo1_str,"(I2)") Histo1
  read(Histo2_str,"(I2)") Histo2
  read(data_str,"(F10.6)") data
  read(DeltaN_str,"(F10.6)") DeltaN
  read(cachestr, "(I3)") cache

  print *, '=========================================================='
  print *, '|| Baur analysis with histogram #', Histo1, '           ||'
  print *, '|| Integrated luminosity :', data,'      ||'
  print *, '|| Theoretical uncertainty :', deltaN,'    ||'
  print *, '|| Delta F_V =', Delf1V,'                                 ||'
  print *, '|| Delta F_VA =', Delf1A,'                                ||'
  print *, '|| Output written to file fort.', cache,'             ||'
  print *, '=========================================================='
  print *, ''
  print *, 'Details:'
  PreFac=PreFac*data

  call GetBinnedDXSec(SM_infile,Histo1,PreFac,NumBins_SM,Valu_SM)
  call GetBinnedDXSec(BSM_infile,Histo1,PreFac,NumBins_BSM,Valu_BSM)
  NumBins=NumBins_SM
  call GetBinnedEvents(Valu_SM,Valu_BSM,NumBins,Events_SM,Events_BSM,NumEventBins)
  if (Histo2 .ne. 0) then      ! then we use two distributions simultaneously
     call GetBinnedDXSec(SM_infile,Histo2,PreFac,NumBins_SM2,Valu_SM2)
     call GetBinnedDXSec(BSM_infile,Histo2,PreFac,NumBins_BSM2,Valu_BSM2)
     NumBins2=NumBins_SM2
     call GetBinnedEvents(Valu_SM2,Valu_BSM2,NumBins2,Events_SM2,Events_BSM2,NumEventBins2)

     do i=1,NumEventBins2
        Events_SM(NumEventBins+i)=Events_SM2(i)
        Events_BSM(NumEventBins+i)=Events_BSM2(i)
     enddo
     NumEventBins=NumEventBins+NumEventBins2

  endif

  print *, 'No. of bins:', NumEventBins
  if (NumEventBins .le. 4) then
     print *, 'NOT EN0UGH EVENTS FOR MEANINGFUL ANALYSIS'
     stop
  endif


  if (trim(operation) .eq. 'GetPVal') then
     call MinChisq(Events_SM,Events_BSM,NumEventBins,DeltaN,chisq,P,sigs)
  elseif (trim(operation) .eq. 'GetLogL') then
     print *, 'Log-likeliood analysis not implemented yet -- to do!'
  endif

   print *, 'chisq = ', chisq
   print *, ''
   print *, 'Output:'
   print *, 'P-value: ', P
   print *, 'sigs',sigs

   write(cache,*) Delf1V,Delf1A,P,sigs


contains 
  
  SUBROUTINE MinChisq(Events_SM,Events_BSM,NumEventBins,DeltaN,chisq,P,sigs)
    implicit none
    integer   :: Events_SM(1:MaxBins),NumEventBins,Events_BSM(1:MaxBins)
    real(8)   :: DeltaN
    real(8)   :: f,fbar_sq,chisq,P,sigs
    integer   :: sum_Events_SM,i,j,k
    
    sum_Events_SM=0
    fbar_sq=0d0

    do i=1,NumEventBins
       print *, i, Events_SM(i),Events_BSM(i)
       sum_Events_SM=sum_Events_SM+Events_SM(i)
       fbar_sq = fbar_sq + 1.d0*Events_BSM(i)**2/Events_SM(i)
    enddo
      
    fbar_sq=fbar_sq/sum_Events_SM
    print *, 'fbar_sq=', fbar_sq

    if ( dsqrt(fbar_sq) .gt. 1d0+DeltaN ) then
       f=1d0+DeltaN
    elseif ( dsqrt(fbar_sq) .lt. 1d0/(1d0+DeltaN) ) then
       f=1d0/(1d0+DeltaN)
    else
       f=dsqrt(fbar_sq)
    endif
    print *, 'f=',f

    chisq=1d0*(NumEventBins-1)
   do j=1,NumEventBins
      chisq=chisq + (Events_BSM(j)*1.0d0 - f*Events_SM(j))**2/(f*Events_SM(j))
   enddo
   
! this should be sorted out already, but put this in here to be sure...
   if (mod(NumEventBins,2).ne.0) then
      print *, 'WARNING: ODD NUMBER OF BINS'
      print *, 'The expression for the P-val is no longer correct'
      stop
   endif
   
 ! now that I have chi-sq, I can get p-value
   P=0d0
   do k=0,NumEventBins/2-1
      P=P+(chisq/2d0)**k/factorial(k)*exp(-chisq/2d0)
   enddo
   sigs=dsqrt(2d0)*inverf(1d0-P)
   print *,'p1=', P
   
   P=pval(chisq/2d0,NumEventBins/2)
   print *,'p2=', P
   pause
   
 end SUBROUTINE MinChisq

! SUBROUTINE GetLogL(Events_SM, Events_BSM, Events_Pred,NumEventBins,DeltaN,L,P,sigs)
!   implicit none
!   integer   :: Events_SM(1:MaxBins),NumEventBins,Events_BSM(1:MaxBins)
!   real(8)   :: DeltaN
!   real(8)   :: L,P,sigs
!   integer   :: SumEvents_SM,SumEvents_BSM
!   integer   :: i
!
!   SumEvents_SM=0
!   SumEvents_BSM=0
!
!   do i=1,NumEventsBins
!      SumEvents_SM =SumEvents_SM +Events_SM(i)
!      SumEvents_BSM=SumEvents_BSM+Events_BSM(i)
!   enddo
!
!   ! now get f
!   ! NB I havent yet verified this, eq (12) of hep-ph/0512262
!   f = 1d0/2d0*(1-deltaN**2*SumEvents_BSM + dsqrt( (1d0-deltaN**2*SumEvents_BSM)**2 + 4d0*deltaN**2*SumEvents_SM ) )
!   
!   if (f .gt. 1d0+deltaN) then
!      f=1d0+deltaN
!   elseif (f .lt. 1d0-deltaN) then
!      f=1d0-deltaN
!   endif
!
!
!   ! this is -2logL
!   lambda = 2d0*
!


  SUBROUTINE GetBinnedEvents(Valu_SM,Valu_BSM,NumBins,Events_SM,Events_BSM,NumEventBins)
    implicit none
    integer,parameter :: MaxBins=200
    real(8)   :: Valu_SM(1:MaxBins),Valu_BSM(1:MaxBins)
    integer   :: NumBins
    real(8)   :: Valu_SM_rebin(1:MaxBins),Valu_BSM_rebin(1:MaxBins)
    integer   :: Events_SM(1:MaxBins),NumEventBins,Events_BSM(1:MaxBins)
    integer   :: i,j,rebin_exp

    j=0
    do i = 1,NumBins
       if (Valu_SM(i) .ge. 5d0 .and. Valu_BSM(i) .ge. 5d0) then
          j=j+1
          Events_SM(j) =floor(Valu_SM(i))
          Events_BSM(j)=floor(Valu_BSM(i))
       endif
    enddo

    rebin_exp=0
    do while ( (j .lt. 5) .and. Numbins/2**rebin_exp .gt. 5) 
       j=0
       rebin_exp=rebin_exp+1
       call REBIN(Valu_SM,Valu_SM_rebin,2**rebin_exp)
       call REBIN(Valu_BSM,Valu_BSM_rebin,2**rebin_exp)
       print *, 'rebinning'
       do i = 1,NumBins/2**rebin_exp
          if (Valu_SM_rebin(i) .ge. 5d0 .and. Valu_BSM_rebin(i) .ge. 5d0) then
             j=j+1
             Events_SM(j) =floor(Valu_SM_rebin(i))
             Events_BSM(j)=floor(Valu_BSM_rebin(i))
          endif
       enddo
    enddo

    if (mod(j,2) .ne. 0) then     ! odd number of bins, set the last to zero
       Events_SM(j) =0
       Events_BSM(j)=0
       j=j-1
    endif

   NumEventBins=j
    
 end SUBROUTINE GetBinnedEvents

  
  SUBROUTINE GetBinnedDXSec(infile,Histo,PreFac,NumBins,Valu)
    implicit none
    character :: infile*(50)
    integer   :: Histo,NumBins 
    integer,parameter :: MaxBins=200
    real(8)   :: Valu(1:MaxBins),PreFac
    integer   :: i 
    character :: dummy*(1)
    integer :: NHisto(1:MaxBins),Hits(1:MaxBins)
    real(8) :: BinSize
    real(8) :: BinVal(1:MaxBins),Error(1:MaxBins)


    open(unit=11,file=trim(infile),form='formatted',access='sequential')

    i=0
    do while(.not.eof(11))  ! loop over all rows

       read(unit=11,fmt="(A)") dummy
       if(dummy(1:1).eq."#") cycle
       backspace(unit=11) ! go to the beginning of the line
      i=i+1
      read(unit=11,fmt=fmt1) NHisto(i),dummy,BinVal(i),dummy,Valu(i),dummy,Error(i),dummy,Hits(i),dummy

      if(NHisto(i).ne.Histo) then
         i=i-1
      else
!         print *,i, BinVal(i)
      endif


   enddo

   NumBins=i
   BinSize=BinVal(2)-BinVal(1)
   Valu  = Valu * PreFac * BinSize

  END SUBROUTINE GetBinnedDXSec


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

  SUBROUTINE REBIN(distr_old,distr_new,div)
    implicit none
    real(8),intent(in) ::  distr_old(:)
    integer, intent(in) :: div
    real(8),intent(out) :: distr_new(1:200)
    integer N,Npr,i,j
! distr with N bins, distr_new has Npr=N/div bins, div >= 1
    
    N = size(distr_old,dim=1)
!    allocate(distr_new(N))
!    if ( mod(N,div) .ne. 0 ) then
!       print *, 'Non integral number of bins!'
!       stop
!    endif

    if ( div .lt. 1 ) then
       print *, 'Final argument should be >= 1'
       stop
    endif

    Npr=N/div
    distr_new = 0d0

    do i = 1,Npr
       do j = 1,div
          distr_new(i) = distr_new(i) + distr_old(div*(i-1)+j)
       enddo
    enddo
    
  end SUBROUTINE REBIN
    
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
       print *, 'INVERF REACHED MAXIMUM 50 NEWTON-RAPHSON ITERATIONS!'
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
    INVERF_KOEFF( 0) =  1.00000000000000D0     
    INVERF_KOEFF( 1) =  1.00000000000000D0     
    INVERF_KOEFF( 2) =  1.16666666666667D0     
    INVERF_KOEFF( 3) =  1.41111111111111D0     
    INVERF_KOEFF( 4) =  1.73373015873016D0     
    INVERF_KOEFF( 5) =  2.14858024691358D0     
    INVERF_KOEFF( 6) =  2.67716623911068D0     
    INVERF_KOEFF( 7) =  3.34814636128128D0     
    INVERF_KOEFF( 8) =  4.19849396342683D0     
    INVERF_KOEFF( 9) =  5.27542268646125D0     
    INVERF_KOEFF(10) =  6.63899509461546D0     
    INVERF_KOEFF(11) =  8.36550450191292D0     
    INVERF_KOEFF(12) =  10.5518020210761D0     
    INVERF_KOEFF(13) =  13.3208081702262D0     
    INVERF_KOEFF(14) =  16.8285215357839D0     
    INVERF_KOEFF(15) =  21.2729253053784D0     
    INVERF_KOEFF(16) =  26.9053027832427D0     
    INVERF_KOEFF(17) =  34.0446123102593D0     
    INVERF_KOEFF(18) =  43.0957486538823D0     
    INVERF_KOEFF(19) =  54.5727422435001D0
    INVERF_KOEFF(20) =  69.1282326361197D0     
    INVERF_KOEFF(21) =  87.5909148253449D0
    INVERF_KOEFF(22) =  111.013117434205D0
          INVERF_KOEFF(23) =    140.731257124519D0     
          INVERF_KOEFF(24) =    178.442657611487D0     
          INVERF_KOEFF(25) =    226.303167593896D0     
          INVERF_KOEFF(26) =    287.051214503983D0     
          INVERF_KOEFF(27) =    364.165459938404D0     
          INVERF_KOEFF(28) =    462.065166575815D0     
          INVERF_KOEFF(29) =    586.364858016595D0     
          INVERF_KOEFF(30) =    744.197995615540D0 
          INVERF_KOEFF(31) =  944.628392281559D0     
          INVERF_KOEFF(32) =  1199.17316418125D0     
          INVERF_KOEFF(33) =  1522.46748209153D0     
          INVERF_KOEFF(34) =  1933.10959970608D0     
          INVERF_KOEFF(35) =  2454.73508332158D0     
          INVERF_KOEFF(36) =  3117.38245243560D0     
          INVERF_KOEFF(37) =  3959.22933516057D0     
          INVERF_KOEFF(38) =  5028.79972696942D0     
          INVERF_KOEFF(39) =  6387.77026383147D0     
          INVERF_KOEFF(40) =  8114.53816822395D0     
          INVERF_KOEFF(41) =  10308.7577173092D0     
          INVERF_KOEFF(42) =  13097.1082841671D0     
          INVERF_KOEFF(43) =  16640.6284810278D0     
          INVERF_KOEFF(44) =   21144.0418418297D0  
!!!    
    if (k .eq. 0) then
       inverf_coeff=1d0
       return
    endif
!    if ( k .le. 44) then
!       INVERF_COEFF = INVERF_KOEFF(k)
!       return
!    endif
    
    do m=0,k-1
       c = c + INVERF_COEFF(m) *  INVERF_COEFF(k-1-m) / ( (m+1)*(2*m+1) )
    enddo
    INVERF_COEFF=c
    return

    
  end function INVERF_COEFF

    
  real function pval(x,n)
    implicit none
! RR Sept 2013
! returns the p-val as defined by P = int^{\infinity}_{x} dx f(x;n)
! with f(x;n) the chisq distribution, f(x;n) = 1/(2^(n/2)) 1/gamma(n/2) x^(n/2-1) exp(-x/2), cf. Cowan eq 2.34
! then P=1/gamma(n/2)*Gamma(n/2,x/2) , with Gamma(n/2,x/2) the upper incomplete gamma function.
! This is given by Gamma(n,x) = gamma(n) - gamma(n,x), where gamma(n,x) is the lower incomplete gamma function
! This is given in series expansion gamma(x,n) = gamma(n) * x^n *exp(-x) * sum_k=0^{\infinity} x^k/gamma(k+n+1) (cf. dlmf.nist.gov/8#PT2)
! Then P(x,n)=1-x^n * exp(-x)* sum_{k=0}^{\infty} x^k/gamma(k+n+1)
! function tested against standard table of P- and chisq values, http://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm

    real(8) ::  x
    integer :: n
    real(8) :: t,old,tol,gamma
    integer :: it

    tol=1d6
    it=0
    pval=0d0
    old=1d6
    
    do while (tol .gt. 1d-10 .and. it .le. 100)
       t=x**it/gamma((n+it+1)*1d0)
       it=it+1
       tol=abs( (t-old)/t )
       old=t
       pval=pval+t
    enddo
    pval=1d0-pval * exp(-x) * x**n
    return

  end function pval
       
    
  

 end program BaurAnalysis

      
   
      

