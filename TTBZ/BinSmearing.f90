


type :: Histogram
    integer :: NBins
    real(8) :: BinSize
    real(8) :: LowVal
    real(8) :: SetScale
    real(8),allocatable :: Value(:)
    real(8),allocatable :: Value2(:)
    integer,allocatable :: Hits(:)
    character :: Info*(50)
    logical :: BinSmearing=.false.
    real(8) :: SmearSigma
end type


 
          Histo(1)%Info    = "pT_lepMinus"
          Histo(1)%NBins   = 40
          Histo(1)%BinSize = 50d0/100d0
          Histo(1)%LowVal  = 0d0
          Histo(1)%SetScale= 100d0
          Histo(1)%BinSmearing= .true.
          Histo(1)%SmearSigma= Histo(1)%BinSize/8d0

          Histo(2)%Info   = "eta_lepMinus"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 0.2d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0
          Histo(2)%BinSmearing= .true.
          Histo(2)%SmearSigma= Histo(2)%BinSize/8d0







SUBROUTINE IntoHisto2(NHisto,BinValue,Value)
use ModMisc
implicit none
integer :: NHisto,NBin,CloserBin
real(8) :: Value,BinValue,ErrorFunct
real(8) :: LowerBinValue,UpperBinValue,CloserBinValue

    NBin = WhichBin(NHisto,BinValue)
    if( IsNaN(Value) ) return
    if( (.not. Histo(NHisto)%BinSmearing) .or. NBin.eq.0 .or. NBin.eq.Histo(NHisto)%NBins+1 ) then
        Histo(NHisto)%Value(NBin) = Histo(NHisto)%Value(NBin)  + Value
        Histo(NHisto)%Value2(NBin)= Histo(NHisto)%Value2(NBin) + Value**2
        Histo(NHisto)%Hits(NBin)  = Histo(NHisto)%Hits(NBin)+1
    else
        LowerBinValue=(NBin-1)*Histo(NHisto)%BinSize + Histo(NHisto)%LowVal
        UpperBinValue=LowerBinValue + Histo(NHisto)%BinSize
        if( BinValue.gt.LowerBinValue+Histo(NHisto)%BinSize/2d0 ) then
           CloserBinValue=UpperBinValue
           CloserBin=NBin+1
           ErrorFunct=erf( dabs(CloserBinValue-BinValue)/Histo(NHisto)%SmearSigma/dsqrt(2d0) )
           Histo(NHisto)%Value(NBin)  = Histo(NHisto)%Value(NBin)  +0.5d0*(1d0+ErrorFunct)*Value
           Histo(NHisto)%Value(NBin+1)= Histo(NHisto)%Value(NBin+1)+0.5d0*(1d0-ErrorFunct)*Value
        else
           CloserBinValue=LowerBinValue
           CloserBin=NBin-1
           ErrorFunct=erf( dabs(CloserBinValue-BinValue)/Histo(NHisto)%SmearSigma/dsqrt(2d0) )
           Histo(NHisto)%Value(NBin)  =Histo(NHisto)%Value(NBin)  +0.5d0*(1d0+ErrorFunct)*Value
           Histo(NHisto)%Value(NBin-1)=Histo(NHisto)%Value(NBin-1)+0.5d0*(1d0-ErrorFunct)*Value
        endif
    endif

RETURN
END SUBROUTINE

