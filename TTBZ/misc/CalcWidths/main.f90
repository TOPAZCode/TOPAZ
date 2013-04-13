PROGRAM IntegrationTemplate
use ModMisc
use ModParameters
use ModIntegrand
use ifport
implicit none
include "vegas_common.f"


   VegasIt1 = 5
   VegasNc1 = 1000000
   integrand = 3


   call coupsm(0)
   call InitVegas()
   call OpenFiles()

   call GetCommandlineArgs()

   call StartVegas()

   call CloseFiles()
END PROGRAM



!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------



SUBROUTINE StartVegas()
use ModParameters
use ModIntegrand
implicit none
include "vegas_common.f"
real(8) :: VG_Result,VG_Error,VG_Chi2
logical :: warmup


if( GridIO.eq.-1 ) then
  readin=.false.
  writeout=.true.
  outgridfile="gridfile"
elseif( GridIO.eq.+1 ) then
  readin=.true.
  writeout=.false.
  ingridfile="gridfile"
else
  readin=.false.
  writeout=.false.
endif
VegasMxDim=mxdim


if( VegasIt0.eq.0 .OR. VegasNc0.eq.0 ) then
   warmup = .false.
   itmx = VegasIt1
   ncall= VegasNc1
else
   itmx = VegasIt0
   ncall= VegasNc0
   warmup = .true.
endif



if(integrand.eq.1) then
    NDim = 2
    call vegas(Integrand_01,VG_Result,VG_Error,VG_Chi2)
elseif(integrand.eq.2) then
    NDim = 2
    call vegas(Integrand_02,VG_Result,VG_Error,VG_Chi2)
elseif(integrand.eq.3) then
    NDim = 5
    call vegas(Integrand_03,VG_Result,VG_Error,VG_Chi2)
elseif(integrand.eq.4) then
    NDim = 5
    call vegas(Integrand_04,VG_Result,VG_Error,VG_Chi2)
else
    print *, "Integrand not available ",integrand
    stop
endif


print *, ""
print *, "Result normalized to 1.46533:",VG_Result/1.46533d0


return
END SUBROUTINE



!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------




SUBROUTINE GetCommandlineArgs()
use ifport
use ModParameters
implicit none
character :: arg*(100)
integer :: NumArgs,NArg


   NumArgs = NArgs()-1
   do NArg=1,NumArgs
    call GetArg(NArg,arg)
    if( arg(1:6).eq."MyArg=" ) then
        read(arg(7:10),*) MyArg
    endif
   enddo



return
END SUBROUTINE




!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------





SUBROUTINE InitVegas()
implicit none
include "vegas_common.f"
integer, parameter :: VegasSeed=19

  idum = -VegasSeed
  xl(1:mxdim) = 0d0
  xu(1:mxdim) = 1d0
  acc = -1d0
  nprn = 1
  readin=.false.
  writeout=.false.

return
END SUBROUTINE






!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------





SUBROUTINE OpenFiles()
use ModParameters
implicit none
character :: filename*(100)

   filename = 'vegas.status'
   open(unit=15,file=trim(filename),form='formatted',access= 'sequential',status='replace')


return
END SUBROUTINE




SUBROUTINE CloseFiles()
implicit none

   close(15)

return
END SUBROUTINE





!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------














