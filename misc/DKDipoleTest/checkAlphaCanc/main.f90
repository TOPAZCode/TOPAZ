PROGRAM IntegrationTemplate
use ModMisc
use ModParameters
use ModIntegrand
use ifport
implicit none
include "vegas_common.f"
integer :: integrand

   m_Top=172d0
   MuRen = m_Top
   alpha_DKfi = 0.01d0

   VegasIt1 = 10
   VegasNc1 = 1000000
   integrand = 1


   call InitVegas()
   call OpenFiles()
   call cpu_time(time_start); print *, "Running",integrand,alpha_DKfi

   call GetCommandlineArgs()

   call StartVegas(integrand)


   call cpu_time(time_end); print *, "Done (",(time_end-time_start)/60d0,") minutes"
   call CloseFiles()
END PROGRAM



!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------



SUBROUTINE StartVegas(integrand)
use ModParameters
use ModIntegrand
implicit none
include "vegas_common.f"
real(8) :: VG_Result,VG_Error,VG_Chi2
logical :: warmup
integer :: integrand


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
    NDim = 8
    call vegas(Integrand_01,VG_Result,VG_Error,VG_Chi2)
elseif(integrand.eq.2) then
    NDim = 5
    call vegas(Integrand_02,VG_Result,VG_Error,VG_Chi2)
else
    print *, "Integrand not available ",integrand
    stop
endif





return
END SUBROUTINE



!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------




SUBROUTINE GetCommandlineArgs(MyArg)
use ifport
use ModParameters
implicit none
character :: arg*(100)
integer :: NumArgs,NArg,MyArg


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














