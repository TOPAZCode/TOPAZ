PROGRAM ttbarjets
use ModParameters
use ModProcess
use ModKinematics
use ModMyRecurrence
use ModJPsiFrag
use ifport
implicit none
include "vegas_common.f"
real(8) :: VG_Result,VG_Error


   call GetCommandlineArgs()
   call Init_cur_2_2f(4)
   call setDim(4,4)
   call InitPDFs()
   call InitParameters()
   call InitHisto()
   call InitPSCuts()
   call InitKirillParameters()
   call InitProcess()
   call InitAmps()
   call InitVegas()
   call InfoPrimAmps("PrimAmpInfo.txt")
   call OpenFiles()
   call WriteParameters(6)   ! stdout
   call WriteParameters(15)  ! vegas status file
   if( TopDecays.eq.5 .or. TopDecays.eq.6 ) call fitFF(MuFrag)

   print *, "Running"
   call cpu_time(time_start)
   call StartVegas(VG_Result,VG_Error)
   call cpu_time(time_end)
!    call WriteHisto(6,it,VG_Result,VG_Error,time_end-time_start)   ! stdout
   call WriteHisto(14,it,VG_Result,VG_Error,time_end-time_start)  ! Histogram file
   call PlotVegas()
   call CloseFiles()
   print *, "Done (",(time_end-time_start)/60d0,") minutes"


END PROGRAM







SUBROUTINE GetCommandlineArgs()
use ModParameters
use ModKinematics
use ifport
implicit none
character :: arg*(100)
character :: env*(31),ColliderStr*(10),ColliderDirStr*(10),ProcessStr*(2),CorrectionStr*(10),FileTag*(50),SeedStr*(6),MuStr*(7),ObsStr*(6)
integer :: NumArgs,NArg,IDipAlpha(1:5),iDKAlpha(1:3),DipAlpha2
logical :: dirresult


   Collider=-1
   Process=-1
   Correction=-1
   ObsSet=-1
   PDFSet=1
   NLOParam=1
   TopDecays=-100
   XTopDecays=-100
   HelSampling=.false.
   m_Top=172d0*GeV
   m_STop=100d0*GeV
   m_Zpr=1500d0*GeV
   Ga_Zpr=m_Zpr*0.01d0
   Q_top = Q_up
   MuRen=m_Top
   MuFac=m_Top
   MuFrag=m_Top
   Fragm_Func_Type=2
   alpha_frag=0.66d0
   beta_frag =12.39d0
   delta_frag=14.97d0
   VegasIt0=-1
   VegasNc0=-1
   VegasIt1=-1
   VegasNc1=-1
   GridIO=0
   Unweighted = .false.
   HistoFile=""
   FileTag=""
   MuStr=""
   ObsStr=""
   GridFile="grid"
   VegasSeed=19
   DKRE_switch=0
   iDipAlpha(1:5)=0
   DipAlpha2=1d0


   NumArgs = NArgs()-1
   do NArg=1,NumArgs
    call GetArg(NArg,arg)
    if( arg(1:9).eq."Collider=" ) then
        read(arg(10:11),*) Collider
    elseif( arg(1:8).eq."Process=" ) then
        read(arg(9:10),*) Process
    elseif( arg(1:11).eq."Correction=" ) then
        read(arg(12:13),*) Correction
    elseif( arg(1:7).eq."ObsSet=" ) then
        read(arg(8:9),*) ObsSet
        write(ObsStr,"(I2)") ObsSet
    elseif( arg(1:5).eq."MTop=" ) then
        read(arg(6:10),*) m_Top
        MuRen=m_Top
        MuFac=m_Top
        MuFrag=m_Top
    elseif( arg(1:6).eq."MStop=" ) then
        read(arg(7:11),*) m_STop
        MuRen=m_STop
        MuFac=m_STop
    elseif( arg(1:5).eq."MZpr=" ) then
        read(arg(6:10),*) m_Zpr
        MuRen=m_Zpr
        MuFac=m_Zpr
    elseif( arg(1:6).eq."GaZpr=" ) then
        read(arg(7:11),*) Ga_Zpr
   elseif( arg(1:8).eq."nGluRad=" ) then
        read(arg(9:10),*) nGluRadContr
    elseif( arg(1:5).eq."XQTop" ) then
        Q_Top = -4d0/3d0
    elseif( arg(1:7).eq."PDFSet=" ) then
        read(arg(8:9),*) PDFSet
    elseif( arg(1:7).eq."FFType=" ) then
        read(arg(8:9),*) Fragm_Func_Type
    elseif( arg(1:5).eq."FFal=" ) then
        read(arg(6:16),*) alpha_frag
    elseif( arg(1:5).eq."FFbe=" ) then
        read(arg(6:16),*) beta_frag
    elseif( arg(1:5).eq."FFde=" ) then
        read(arg(6:16),*) delta_frag
    elseif( arg(1:9).eq."NLOParam=" ) then
        read(arg(10:11),*) NLOParam
    elseif( arg(1:6).eq."MuRen=" ) then
        read(arg(7:10),*) MuRen
!         if( MuRen.ne.m_Top ) write(MuStr,"(F4.2)") MuRen
    elseif( arg(1:6).eq."MuFac=" ) then
        read(arg(7:10),*) MuFac
    elseif( arg(1:7).eq."MuFrag=" ) then
        read(arg(8:11),*) MuFrag
    elseif( arg(1:6).eq."TopDK=" ) then
        read(arg(7:9),*) TopDecays
    elseif( arg(1:7).eq."XTopDK=" ) then
        read(arg(8:10),*) XTopDecays
    elseif( arg(1:9).eq."VegasIt0=" ) then
        read(arg(10:11),*) VegasIt0
    elseif( arg(1:9).eq."VegasIt1=" ) then
        read(arg(10:11),*) VegasIt1
    elseif( arg(1:9).eq."VegasNc0=" ) then
        read(arg(10:20),*) VegasNc0
    elseif( arg(1:9).eq."VegasNc1=" ) then
        read(arg(10:20),*) VegasNc1
    elseif( arg(1:10).eq."VegasSeed=" ) then
        read(arg(11:17),*) VegasSeed
    elseif( arg(1:9).eq."GridFile=" ) then
        read(arg(10:41),*) GridFile
    elseif( arg(1:7).eq."GridIO=" ) then
        read(arg(8:10),*) GridIO
    elseif( arg(1:10).eq."HistoFile=" ) then
        read(arg(11:41),*) HistoFile
    elseif( arg(1:8).eq."FileTag=" ) then
        read(arg(9:59),*)  FileTag
        if( FileTag.eq."." ) FileTag=""
    elseif( arg(1:9).eq."DipAlpha=" ) then
        read(arg(10:10),*) iDipAlpha(1)
        read(arg(11:11),*) iDipAlpha(2)
        read(arg(12:12),*) iDipAlpha(3)
        read(arg(13:13),*) iDipAlpha(4)
        read(arg(14:14),*) iDipAlpha(5)
        alpha_ii = 10d0**(-iDipAlpha(1))
        alpha_if = 10d0**(-iDipAlpha(2))
        alpha_fi = 10d0**(-iDipAlpha(3))
        alpha_ff = 10d0**(-iDipAlpha(4))
        alpha_DK = 10d0**(-iDipAlpha(5))
    elseif( arg(1:8).eq."DKAlpha=" ) then
        read(arg(9:9),*)   iDKAlpha(1)
        read(arg(10:10),*) iDKAlpha(2)
        read(arg(11:11),*) iDKAlpha(3)
    elseif( arg(1:10).eq."DipAlpha2=" ) then
        read(arg(11:11),*) DipAlpha2
    elseif( arg(1:8).eq."HelSamp=" ) then
        read(arg(9:9),*) HelSampling
    elseif( arg(1:5).eq."DKRE=" ) then
        read(arg(6:7),*) DKRE_switch
    endif
   enddo
   write(MuStr,"(F5.2)") MuRen

   if( DipAlpha2.eq.0d0 ) then
       print *, "DipAlpha2 cannot ne zero"
       !stop
   endif
   if (alpha_ii.ne.1d0) alpha_ii = DipAlpha2 * alpha_ii
   if (alpha_if.ne.1d0) alpha_if = DipAlpha2 * alpha_if
   if (alpha_fi.ne.1d0) alpha_fi = DipAlpha2 * alpha_fi
   if (alpha_ff.ne.1d0) alpha_ff = DipAlpha2 * alpha_ff


   alpha_DKTfi = alpha_DK; alpha_DKTff = alpha_DK; alpha_DKWff = alpha_DK;
   if(iDKAlpha(1).ne.0) then
      alpha_DKTfi = dble(iDKAlpha(1)) * alpha_DKTfi
   endif
   if(iDKAlpha(2).ne.0) then
      alpha_DKTff = dble(iDKAlpha(2)) * alpha_DKTff
   endif
   if(iDKAlpha(3).ne.0) then
      alpha_DKWff = dble(iDKAlpha(3)) * alpha_DKWff
   endif


    if(Collider.eq.-1 .or. Process.eq.-1 .or. Correction.eq.-1 .or. TopDecays.eq.-100 .or. ObsSet.eq.-1) then
          print *, "not enough input parameter"
          print *, "required: Collider,Process,Correction,TopDK,ObsSet"
          print *, "Collider:    1=LHC(14TeV), 11=LHC(7TeV), 12=LHC(8TeV), 13=LHC(13TeV), 2=Tevatron"
          print *, "Process:     see ProcessInfo.txt"
          print *, "Correction:  0=LO, 1=VI, 2=RE, 3=ID, 4=VI in top decay, 5=RE in top decay"
          print *, "TopDK:       0=stable tops, 1=dilept, 2=full hadr, 3=l-&jets, 4=l+&jets"
          print *, "ObsSet:      see mod_Kinematics.f90"
          stop
    endif
    if( Process.ge.41 .and. Process.le.59 ) then
          if( XTopDecays.eq.-100 ) then 
              print *, "not enough input parameter"
              print *, "required: XTopDK:      0=stable, 1=Htop-->vector+top, 2=HTop-->scalar+top, 3=Stop-->Chi0+top"
          endif
    endif

    if((Correction.eq.3 .or. TopDecays.eq.5 .or. TopDecays.eq.6) .and. HelSampling) then
          print *, "Helicity sampling is not allowed here. This is simpliy due to assignments of random numbers in mod_CrossSection.f90"
    endif

    if(Process.lt.10) then
      write(ProcessStr,"(I1)") Process
      ProcessStr="0"//trim(ProcessStr)
    else
      write(ProcessStr,"(I2)") Process
    endif

    if( Collider.eq.1 ) then
         AlgoType = -1       ! anti-kT
    elseif( Collider.eq.2 ) then
         AlgoType = +1       ! kT
    endif


    if(Correction.eq.0) then
        if(NLOParam.le.1) then
          CorrectionStr = "LOLO"
        elseif(NLOParam.eq.2) then
          CorrectionStr = "LO"
        endif
    elseif(Correction.eq.1) then
        CorrectionStr = "1L"
    elseif(Correction.eq.2) then
        CorrectionStr = "RE"
    elseif(Correction.eq.3) then
        CorrectionStr = "ID"
    elseif(Correction.eq.4) then
        CorrectionStr = "DK1L"
    elseif(Correction.eq.5 .and. DKRE_switch.eq.1) then
        CorrectionStr = "DKRE_a"
    elseif(Correction.eq.5 .and. DKRE_switch.eq.2) then
        CorrectionStr = "DKRE_b"
    elseif(Correction.eq.5 .and. DKRE_switch.eq.0) then
        CorrectionStr = "DKRE"
    endif


    if(VegasSeed.eq.19) then
      SeedStr=""
    elseif(VegasSeed.lt.10) then
      write(SeedStr,"(I1)") VegasSeed
      SeedStr="_00"//trim(SeedStr)
    elseif(VegasSeed.lt.100) then
      write(SeedStr,"(I2)") VegasSeed
      SeedStr="_0"//trim(SeedStr)
    elseif(VegasSeed.lt.1000) then
      write(SeedStr,"(I3)") VegasSeed
      SeedStr="_"//trim(SeedStr)
    else
      print *, "Vegas seed too big"
      stop
    endif

   if(ObsSet.lt.10) then
      write(ObsStr,"(I1)") ObsSet
      ObsStr="0"//trim(ObsStr)
   endif

    if( Collider.eq.1 ) then
        ColliderDirStr="LHC14"
        ColliderStr="LHC"
    elseif( Collider.eq.11 ) then
        ColliderDirStr="LHC7"
        ColliderStr="LHC"
    elseif( Collider.eq.12 ) then
        ColliderDirStr="LHC8"
        ColliderStr="LHC"
    elseif( Collider.eq.13 ) then
        ColliderDirStr="LHC13"
        ColliderStr="LHC"
    elseif( Collider.eq.2 ) then
        ColliderDirStr="TEV"
        ColliderStr="TEV"
    endif


    if( Q_top.ne.Q_up ) then
        MuStr=trim(MuStr)//"_XQ"
    endif

    dirresult = makedirqq("./"//trim(ColliderDirStr)//"_"//trim(ObsStr)//"_"//trim(MuStr))
    dirresult = makedirqq("./"//trim(ColliderDirStr)//"_"//trim(ObsStr)//"_"//trim(MuStr)//"/"//trim(ProcessStr))
    if(dirresult) print *, "created directory "//"./"//trim(ColliderDirStr)//"_"//trim(ObsStr)//"_"//trim(MuStr)//"/"//trim(ProcessStr)


    if( ObsSet.eq.8 ) then!   spin correlations with R
            if(TopDecays.eq.+1) FileTag=trim(FileTag)//"_c"
            if(TopDecays.eq.-1) FileTag=trim(FileTag)//"_u"
            print *, "Remember: MRST PDFs, Ellis-Soper Recombination"
    endif


    if(HistoFile.eq."") then
        HistoFile = "./"//trim(ColliderDirStr)//"_"//trim(ObsStr)//"_"//trim(MuStr)//"/"//trim(ProcessStr)//"/"//trim(ColliderStr)//"."//trim(ProcessStr)//"."//trim(CorrectionStr)//trim(SeedStr)//trim(FileTag)
    endif

   GridFile = "./"//trim(ColliderDirStr)//"_"//trim(ObsStr)//"_"//trim(MuStr)//"/"//trim(ProcessStr)//"/"//trim(ColliderStr)//"."//trim(ProcessStr)//".1L."//trim(GridFile)

return
END SUBROUTINE






SUBROUTINE WriteParameters(TheUnit)
use ModParameters
use ModKinematics
implicit none
integer TheUnit


!DEC$ IF(_CheckMomenta .EQ.1)
    print *, "_CheckMomenta activated!"
!DEC$ ENDIF


   write(TheUnit,"(A)") "# Program parameters:"
   write(TheUnit,"(A,I2,A,F8.3,A)") "# Collider=",Collider," (",Collider_Energy*1d-1," TeV)"
   write(TheUnit,"(A,I2)") "# Process=",Process
   write(TheUnit,"(A,I2)") "# Master Process=",MasterProcess
   write(TheUnit,"(A,I2)") "# Correction=",Correction
   write(TheUnit,"(A,I2,A)") "# PDF Set=",PDFSet,PDFSetString
   write(TheUnit,"(A,I2)") "# NLO Parameter=",NLOParam
   write(TheUnit,"(A,I2)") "# Top Decays=",TopDecays
   write(TheUnit,"(A,F7.2,A)") "# MuRen=",MuRen*100," GeV"
   write(TheUnit,"(A,F7.2,A)") "# MuFac=",MuFac*100," GeV"
   if( Q_top.ne.Q_up ) write(TheUnit,"(A,F13.9)") "# Q_top=",Q_top
   if( Correction.eq.2 .or. Correction.eq.3 .or. Correction.eq.4 .or. Correction.eq.5 ) then
        write(TheUnit,"(A,PE9.2,PE9.2,PE9.2,PE9.2,PE9.2,A)") "# Alpha Parameters (ii if fi ff DK)=(",alpha_ii,alpha_if,alpha_fi,alpha_ff,alpha_DK,")"
        write(TheUnit,"(A,PE9.2,PE9.2,PE9.2,A)") "# DK Alpha Parameters (Tfi Tff Wff)=(",alpha_DKTfi,alpha_DKTff,alpha_DKWff,")"
        write(TheUnit,"(A,F13.6)") "# Kappa Parameter  kappa_ff=",kappa_ff
   endif
   if( TOPDECAYS.EQ.5 .OR. TOPDECAYS.EQ.6  ) then
      write(TheUnit,"(A,F7.2,A)") "# MuFrag=",MuFrag*100," GeV"
      write(TheUnit,"(A,I3)") "# Fragmentation function type=",Fragm_Func_Type
      write(TheUnit,"(A,F10.5,F10.5,F10.5)") "# parameters (alpha,beta,delta)=",alpha_frag,beta_frag,delta_frag
   endif

    if(MuRen.ne.MuFac) then
       write(TheUnit,*) "# MuRen.ne.MuFac: check that this is correctly implemented in the dipole routines!"
    endif

    write(TheUnit,"(A,F13.9,A)") "# alpha_s(MuRen)=",alpha_s*RunAlphaS(NLOParam,MuRen)
    write(TheUnit,"(A,F13.9,A)") "# 1/alpha_em=",1d0/alpha
    if(NLOParam.eq.1) then
       write(TheUnit,"(A)") "# one loop running"
    elseif(NLOParam.eq.2) then
       write(TheUnit,"(A)") "# two loop running"
    else
       write(TheUnit,*) "# no alpha_s running"
    endif
    if( ObsSet.ge.60 .and. ObsSet.le.69 ) then
        write(TheUnit,'(A,F10.5)') "# m(Zpr)=",m_Zpr*100d0
        write(TheUnit,'(A,F10.5)') "# Gamma(Zpr)=",Ga_Zpr*100d0
    endif
    write(TheUnit,'(A,F8.3)') "# m(top)=",m_Top *100d0
    write(TheUnit,"(A,F10.6)") "# Width expansion factor=",WidthExpansion
    write(TheUnit,"(A,F10.6,A)") "# Gamma_Top(LO) =",Ga_Top(0)*100d0," GeV"
    write(TheUnit,"(A,F10.6,A)") "# Gamma_Top(NLO)=",(Ga_Top(0)+Ga_Top(1))*100d0," GeV"
    write(TheUnit,"(A,F10.6,A)") "# Gamma_W(LO) =",Ga_W(0)*100d0," GeV"
    write(TheUnit,"(A,F10.6,A)") "# Gamma_W(NLO)=",(Ga_W(0)+Ga_W(1))*100d0," GeV"
    if( AlgoType.eq.1 ) then
        write(TheUnit,"(A,A)") "# Jet Algorithm= kT"
    elseif( AlgoType.eq.-1 ) then
        write(TheUnit,"(A,A)") "# Jet Algorithm= anti kT"
    else
        write(TheUnit,"(A)") "# Jet Algorithm= UNKNOWN"
    endif
    if( RecombPrescr.eq.0 ) then
        write(TheUnit,"(A,A)") "# Recombination scheme= 4-vector addition"
    elseif( RecombPrescr.eq.1 ) then
        write(TheUnit,"(A,A)") "# Recombination scheme= Ellis-Soper"
    else
        write(TheUnit,"(A)") "# Recombination scheme= UNKNOWN"
    endif

    if(HelSampling) then
       write(TheUnit,"(A)") "# helicity sampling"
    else
       write(TheUnit,"(A)") "# no helicity sampling"
    endif


   write(TheUnit,'(A,F8.3)') "# cuts:"
   write(TheUnit,'(A,F8.3)') "# pT_jet_cut= ",pT_jet_cut*100d0
   write(TheUnit,'(A,F8.3)') "# pT_hardestjet_cut= ",pT_hardestjet_cut*100d0
   write(TheUnit,'(A,F8.3)') "# eta_jet_cut= ",eta_jet_cut
   write(TheUnit,'(A,F8.3)') "# eta_sepa_cut= ",eta_sepa_cut
   write(TheUnit,'(A,F8.3)') "# pT_bjet_cut= ",pT_bjet_cut*100d0
   write(TheUnit,'(A,F8.3)') "# eta_bjet_cut= ",eta_bjet_cut
   write(TheUnit,'(A,F8.3)') "# Rsep_jet= ",Rsep_jet
   write(TheUnit,'(A,F8.3)') "# pT_lep_cut= ",pT_lep_cut*100d0
   write(TheUnit,'(A,F8.3)') "# eta_lep_cut= ",eta_lep_cut
   write(TheUnit,'(A,F8.3)') "# pT_miss_cut= ",pT_miss_cut*100d0
   write(TheUnit,'(A,F8.3)') "# MInv_jets_cut= ",MInv_jets_cut*100d0
   write(TheUnit,'(A,F8.3)') "# HT_cut= ",HT_cut*100d0
   write(TheUnit,'(A,F8.3)') "# pt_photon_cut= ",pt_pho_cut*100d0
   write(TheUnit,'(A,F8.3)') "# eta_photon_cut= ",eta_pho_cut
   write(TheUnit,'(A,F8.3)') "# Rsep_photonjet= ",Rsep_Pj
   write(TheUnit,'(A,F8.3)') "# Rsep_photonbjet= ",Rsep_Pbj
   write(TheUnit,'(A,F8.3)') "# Rsep_photonlepton= ",Rsep_Plep

   write(TheUnit,'(A)') "#"

   write(TheUnit,"(A,I15)") "# VegasSeed=",VegasSeed
   write(TheUnit,"(A,I15)") "# VegasIt0 =",VegasIt0
   write(TheUnit,"(A,I15)") "# VegasNc0 =",VegasNc0
   write(TheUnit,"(A,I15)") "# VegasIt1 =",VegasIt1
   write(TheUnit,"(A,I15)") "# VegasNc1 =",VegasNc1
   write(TheUnit,"(A)") "# Histo file:        "//trim(HistoFile)//'.dat'
   write(TheUnit,"(A)") "# Vegas status file: "//trim(HistoFile)//'.status'
   write(TheUnit,"(A)") "# Histo status file: "//trim(HistoFile)//'.tmp_histo'


END SUBROUTINE





!DEC$ IF(_UseCuba .EQ.0)
SUBROUTINE StartVegas(VG_Result,VG_Error)
use ModMisc
use ModCrossSection_TTB
use ModCrossSection_TTBJ
use ModCrossSection_TTBP
use ModCrossSection_TTBETmiss
use ModCrossSection_ZprimeTTB
use ModKinematics
use ModParameters
implicit none
include "vegas_common.f"
real(8) :: VG_Result,VG_Error,VG_Chi2
logical :: warmup


if( GridIO.eq.-1 ) then
  readin=.false.
  writeout=.true.
  outgridfile=GridFile(1:72)
elseif( GridIO.eq.+1 ) then
  readin=.true.
  writeout=.false.
  ingridfile=GridFile(1:72)
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

IF( MASTERPROCESS.EQ.0 ) THEN
IF( CORRECTION   .EQ.1 ) THEN
  call vegas(EvalCS_1L_gggggg,VG_Result,VG_Error,VG_Chi2)
ENDIF
ENDIF


IF( MASTERPROCESS.EQ.1 ) THEN
IF( CORRECTION.LE.1 .AND. PROCESS.EQ.1 .AND. TOPDECAYS.NE.101) THEN
  call vegas(EvalCS_1L_ttbgg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_1L_ttbgg,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.LE.1 .AND. PROCESS.EQ.1 .AND. TOPDECAYS.EQ.101) THEN! this is experimental: ttbar off-shell production at LO
  ndim=14
  call vegas(EvalCS_LO_bbbWWgg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_LO_bbbWWgg,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.LE.1 .AND. PROCESS.EQ.21 ) THEN
  call vegas(EvalCS_DKP_1L_ttbgg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_DKP_1L_ttbgg,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.LE.1 .AND. PROCESS.EQ.33 ) THEN
  call vegas(EvalCS_DKJ_1L_ttbgg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_DKJ_1L_ttbgg,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.LE.1 .AND. PROCESS.EQ.41 ) THEN
  m_Top = m_SMTop!   restore correct top mass value since initialization is done (see InitProcess(PROCESS.EQ.41))
  call vegas(EvalCS_1L_HtHtbgg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_1L_HtHtbgg,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.3 .AND. PROCESS.EQ.5 ) THEN
  call vegas(EvalCS_1L_ttbgg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_1L_ttbgg,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.3 .AND. PROCESS.EQ.29 ) THEN
  call vegas(EvalCS_DKP_1L_ttbgg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_DKP_1L_ttbgg,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.3 .AND. PROCESS.EQ.37 ) THEN
  call vegas(EvalCS_DKJ_1L_ttbgg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_DKJ_1L_ttbgg,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.4 .AND. PROCESS.EQ.1 ) THEN
  if( TOPDECAYS.GT.0 ) call vegas(EvalCS_NLODK_ttb,     VG_Result,VG_Error,VG_Chi2)
  if( TOPDECAYS.LT.0 ) call vegas(EvalCS_NLODK_ttb_noSC,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   if( TOPDECAYS.GT.0 ) call vegas1(EvalCS_NLODK_ttb,     VG_Result,VG_Error,VG_Chi2)
   if( TOPDECAYS.LT.0 ) call vegas1(EvalCS_NLODK_ttb_noSC,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.4 .AND. PROCESS.EQ.21 ) THEN
  call vegas(EvalCS_NLODKP_ttb,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_NLODKP_ttb,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.4 .AND. PROCESS.EQ.33 ) THEN
  call vegas(EvalCS_1LDKJ_ttb,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_1LDKJ_ttb,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.5 .AND. PROCESS.EQ.1 ) THEN
  if( TOPDECAYS.GT.0 ) call vegas(EvalCS_NLODK_ttb,     VG_Result,VG_Error,VG_Chi2)
  if( TOPDECAYS.LT.0 ) call vegas(EvalCS_NLODK_ttb_noSC,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   if( TOPDECAYS.GT.0 ) call vegas1(EvalCS_NLODK_ttb,     VG_Result,VG_Error,VG_Chi2)
   if( TOPDECAYS.LT.0 ) call vegas1(EvalCS_NLODK_ttb_noSC,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.5 .AND. PROCESS.EQ.21 ) THEN
  call vegas(EvalCS_NLODKP_ttb,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_NLODKP_ttb,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.5 .AND. PROCESS.EQ.33 ) THEN
  call vegas(EvalCS_REDKJ_ttb,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_REDKJ_ttb,VG_Result,VG_Error,VG_Chi2)
  endif
ENDIF
ENDIF




IF( MASTERPROCESS.EQ.2 ) THEN
IF( CORRECTION.LE.1 .AND. PROCESS.EQ.2 ) THEN
        call vegas(EvalCS_1L_ttbqqb,VG_Result,VG_Error,VG_Chi2)
        if( warmup ) then
        itmx = VegasIt1
        ncall= VegasNc1
        call InitHisto()
        call vegas1(EvalCS_1L_ttbqqb,VG_Result,VG_Error,VG_Chi2)
        endif
ELSEIF( CORRECTION.LE.1 .AND. PROCESS.EQ.23 ) THEN
        call vegas(EvalCS_DKP_1L_ttbqqb,VG_Result,VG_Error,VG_Chi2)
        if( warmup ) then
        itmx = VegasIt1
        ncall= VegasNc1
        call InitHisto()
        call vegas1(EvalCS_DKP_1L_ttbqqb,VG_Result,VG_Error,VG_Chi2)
        endif
ELSEIF( CORRECTION.LE.1 .AND. PROCESS.EQ.34 ) THEN
        call vegas(EvalCS_DKJ_1L_ttbqqb,VG_Result,VG_Error,VG_Chi2)
        if( warmup ) then
        itmx = VegasIt1
        ncall= VegasNc1
        call InitHisto()
        call vegas1(EvalCS_DKJ_1L_ttbqqb,VG_Result,VG_Error,VG_Chi2)
        endif
ELSEIF( CORRECTION.LE.1 .AND. PROCESS.EQ.42 ) THEN
  m_Top = m_SMTop!  restore correct top mass value since initialization is done (see InitProcess(PROCESS.EQ.42))
  call vegas(EvalCS_1L_HtHtbqqb,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_1L_HtHtbqqb,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.3 .AND. PROCESS.LE.6 ) THEN
  call vegas(EvalCS_1L_ttbqqb,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_1L_ttbqqb,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.3 .AND. (PROCESS.EQ.25.OR.PROCESS.EQ.27.OR.PROCESS.EQ.31) ) THEN
  call vegas(EvalCS_DKP_1L_ttbqqb,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_DKP_1L_ttbqqb,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.3 .AND. (PROCESS.EQ.35.OR.PROCESS.EQ.36.OR.PROCESS.EQ.38) ) THEN
  call vegas(EvalCS_DKJ_1L_ttbqqb,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_DKJ_1L_ttbqqb,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.4 .AND. PROCESS.EQ.2 ) THEN
  if( TOPDECAYS.GT.0 ) call vegas(EvalCS_NLODK_ttb,     VG_Result,VG_Error,VG_Chi2)
  if( TOPDECAYS.LT.0 ) call vegas(EvalCS_NLODK_ttb_noSC,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   if( TOPDECAYS.GT.0 ) call vegas1(EvalCS_NLODK_ttb,     VG_Result,VG_Error,VG_Chi2)
   if( TOPDECAYS.LT.0 ) call vegas1(EvalCS_NLODK_ttb_noSC,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.4 .AND. PROCESS.EQ.23 ) THEN
  call vegas(EvalCS_NLODKP_ttb,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_NLODKP_ttb,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.4 .AND. PROCESS.EQ.34 ) THEN
  call vegas(EvalCS_1LDKJ_ttb,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_1LDKJ_ttb,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.5 .AND. PROCESS.EQ.2 ) THEN
  if( TOPDECAYS.GT.0 ) call vegas(EvalCS_NLODK_ttb,     VG_Result,VG_Error,VG_Chi2)
  if( TOPDECAYS.LT.0 ) call vegas(EvalCS_NLODK_ttb_noSC,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   if( TOPDECAYS.GT.0 ) call vegas1(EvalCS_NLODK_ttb,     VG_Result,VG_Error,VG_Chi2)
   if( TOPDECAYS.LT.0 ) call vegas1(EvalCS_NLODK_ttb_noSC,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.5 .AND. PROCESS.EQ.23 ) THEN
  call vegas(EvalCS_NLODKP_ttb,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_NLODKP_ttb,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.5 .AND. PROCESS.EQ.34 ) THEN
  call vegas(EvalCS_REDKJ_ttb,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_REDKJ_ttb,VG_Result,VG_Error,VG_Chi2)
  endif
ENDIF
ENDIF





IF( MASTERPROCESS.EQ.3 ) THEN
IF( CORRECTION   .LE.1 ) THEN
  call vegas(EvalCS_1L_ttbggg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_1L_ttbggg,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.2 .AND. PROCESS.EQ.5) THEN
  call vegas(EvalCS_Real_ttbggg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_Real_ttbggg,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.2 .AND. PROCESS.EQ.29) THEN
  call vegas(EvalCS_DKP_Real_ttbggg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_DKP_Real_ttbggg,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.2 .AND. PROCESS.EQ.37) THEN
  call vegas(EvalCS_DKJ_Real_ttbggg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_DKJ_Real_ttbggg,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION .EQ.3 ) THEN
  call vegas(EvalCS_1L_ttbggg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_1L_ttbggg,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION .EQ.4 ) THEN
  call vegas(EvalCS_1LDK_ttbggg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_1LDK_ttbggg,VG_Result,VG_Error,VG_Chi2)
  endif
ENDIF
ENDIF




IF( MASTERPROCESS.EQ.4 ) THEN
IF( CORRECTION.EQ.0 ) THEN
  call vegas(EvalCS_1L_ttbqqbg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_1L_ttbqqbg,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.2 .AND. (PROCESS.EQ.3 .OR. PROCESS.EQ.4 .OR. PROCESS.EQ.6) ) THEN
  call vegas(EvalCS_Real_ttbqqbg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_Real_ttbqqbg,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.2 .AND. (PROCESS.EQ.25 .OR. PROCESS.EQ.27 .OR. PROCESS.EQ.31) ) THEN
  call vegas(EvalCS_DKP_Real_ttbqqbg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_DKP_Real_ttbqqbg,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.2 .AND. (PROCESS.EQ.35 .OR. PROCESS.EQ.36 .OR. PROCESS.EQ.38) ) THEN
  call vegas(EvalCS_DKJ_Real_ttbqqbg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_DKJ_Real_ttbqqbg,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.1 ) THEN
  call vegas(EvalCS_1L_ttbqqbg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_1L_ttbqqbg,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.3 ) THEN
  call vegas(EvalCS_1L_ttbqqbg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_1L_ttbqqbg,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.4 ) THEN
  call vegas(EvalCS_1LDK_ttbqqbg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_1LDK_ttbqqbg,VG_Result,VG_Error,VG_Chi2)
  endif
ENDIF
ENDIF



IF( MASTERPROCESS.EQ.5 ) THEN
IF( CORRECTION   .EQ.2 ) THEN
  call vegas(EvalCS_Real_ttbgggg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_Real_ttbgggg,VG_Result,VG_Error,VG_Chi2)
  endif
ENDIF
ENDIF



IF( MASTERPROCESS.EQ.6 ) THEN
IF( CORRECTION   .EQ.2 ) THEN
  call vegas(EvalCS_Real_ttbqqbgg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_Real_ttbqqbgg,VG_Result,VG_Error,VG_Chi2)
  endif
ENDIF
ENDIF



IF( MASTERPROCESS.EQ.7 ) THEN
IF( CORRECTION   .EQ.2 ) THEN
  call vegas(EvalCS_Real_ttbqqbqqb,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_Real_ttbqqbqqb,VG_Result,VG_Error,VG_Chi2)
  endif
ENDIF
ENDIF




IF( MASTERPROCESS.EQ.8 ) THEN
IF( CORRECTION   .LE.1 ) THEN
  call vegas(EvalCS_1L_ttbggp,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_1L_ttbggp,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.3 ) THEN
  call vegas(EvalCS_1L_ttbggp,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_1L_ttbggp,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.4 ) THEN
  call vegas(EvalCS_NLODK_ttbp,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_NLODK_ttbp,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.5 ) THEN
  call vegas(EvalCS_NLODK_ttbp,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_NLODK_ttbp,VG_Result,VG_Error,VG_Chi2)
  endif
ELSE
  call Error("this correction is not available")
ENDIF
ENDIF


IF( MASTERPROCESS.EQ.9 ) THEN
IF( CORRECTION   .LE.1 ) THEN
  call vegas(EvalCS_1L_ttbqqbp,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_1L_ttbqqbp,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.3 ) THEN
  call vegas(EvalCS_1L_ttbqqbp,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_1L_ttbqqbp,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.4 ) THEN
  call vegas(EvalCS_NLODK_ttbp,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_NLODK_ttbp,VG_Result,VG_Error,VG_Chi2)
  endif
ELSEIF( CORRECTION.EQ.5 ) THEN
  call vegas(EvalCS_NLODK_ttbp,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_NLODK_ttbp,VG_Result,VG_Error,VG_Chi2)
  endif
ENDIF
ENDIF



IF( MASTERPROCESS.EQ.10 ) THEN
IF( CORRECTION.EQ.2 ) THEN
  call vegas(EvalCS_Real_ttbgggp,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_Real_ttbgggp,VG_Result,VG_Error,VG_Chi2)
  endif
ENDIF
ENDIF




IF( MASTERPROCESS.EQ.11 ) THEN
IF( CORRECTION.EQ.2 ) THEN
  call vegas(EvalCS_Real_ttbqqbgp,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_Real_ttbqqbgp,VG_Result,VG_Error,VG_Chi2)
  endif
ENDIF
ENDIF




IF( MASTERPROCESS.EQ.12 ) THEN
IF( CORRECTION.LE.1 .AND. PROCESS.EQ.51 ) THEN
  call vegas(EvalCS_1L_ststbgg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
    itmx = VegasIt1
    ncall= VegasNc1
    call InitHisto()
    call vegas1(EvalCS_1L_ststbgg,VG_Result,VG_Error,VG_Chi2)
  endif

ELSEIF( CORRECTION.EQ.3 .AND. PROCESS.EQ.55 ) THEN
  call vegas(EvalCS_1L_ststbgg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
    itmx = VegasIt1
    ncall= VegasNc1
    call InitHisto()
    call vegas1(EvalCS_1L_ststbgg,VG_Result,VG_Error,VG_Chi2)
  endif
ENDIF
ENDIF




IF( MASTERPROCESS.EQ.13 ) THEN
IF( CORRECTION.LE.1 .AND. PROCESS.EQ.52 ) THEN
  call vegas(EvalCS_1L_ststbqqb,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
    itmx = VegasIt1
    ncall= VegasNc1
    call InitHisto()
    call vegas1(EvalCS_1L_ststbqqb,VG_Result,VG_Error,VG_Chi2)
  endif

ELSEIF( CORRECTION.EQ.3 .AND. (PROCESS.EQ.56 .OR. PROCESS.EQ.53 .OR. PROCESS.EQ.54 ) ) THEN
  call vegas(EvalCS_1L_ststbqqb,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
    itmx = VegasIt1
    ncall= VegasNc1
    call InitHisto()
    call vegas1(EvalCS_1L_ststbqqb,VG_Result,VG_Error,VG_Chi2)
  endif
ENDIF
ENDIF




IF( MASTERPROCESS.EQ.14 ) THEN
IF( CORRECTION.EQ.2 .AND. PROCESS.EQ.55 ) THEN
  call vegas(EvalCS_Real_ststbggg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
    itmx = VegasIt1
    ncall= VegasNc1
    call InitHisto()
    call vegas1(EvalCS_Real_ststbggg,VG_Result,VG_Error,VG_Chi2)
  endif
ENDIF
ENDIF





IF( MASTERPROCESS.EQ.15 ) THEN
IF( CORRECTION.EQ.2 .AND. (PROCESS.EQ.53 .OR. PROCESS.EQ.54 .OR. PROCESS.EQ.56 ) ) THEN
  call vegas(EvalCS_Real_ststbqqbg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
    itmx = VegasIt1
    ncall= VegasNc1
    call InitHisto()
    call vegas1(EvalCS_Real_ststbqqbg,VG_Result,VG_Error,VG_Chi2)
  endif
ENDIF
ENDIF








IF( MASTERPROCESS.EQ.16 ) THEN
IF( CORRECTION.LE.1 .AND. (PROCESS.EQ.56 .OR. PROCESS.EQ.59) ) THEN
  call vegas(EvalCS_1L_ststbgggg,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
    itmx = VegasIt1
    ncall= VegasNc1
    call InitHisto()
    call vegas1(EvalCS_1L_ststbgggg,VG_Result,VG_Error,VG_Chi2)
  endif
ENDIF
ENDIF



!!! Zprime section !!!

IF( MASTERPROCESS.EQ.62 ) THEN
IF( (CORRECTION.LE.1) .AND. PROCESS.EQ.62 ) THEN
   print *, 'LO / Virt Zprime'
   call qlinit()
   call vegas(EvalCS_1L_Zprime_ttbqqb,VG_Result,VG_Error,VG_Chi2)
   if( warmup ) then
      itmx = VegasIt1
      ncall= VegasNc1
      call InitHisto()
      call vegas1(EvalCS_1L_Zprime_ttbqqb,VG_Result,VG_Error,VG_Chi2)
   endif
ELSEIF ( CORRECTION.EQ.3 .AND. (PROCESS.EQ.63 .OR. PROCESS.EQ.64 .OR. PROCESS.EQ.66)) THEN
   print *, 'Zprime Int Dipoles'
   call qlinit()
   call vegas(EvalCS_1L_Zprime_ttbqqb,VG_Result,VG_Error,VG_Chi2)
   if( warmup ) then
      itmx = VegasIt1
      ncall= VegasNc1
      call InitHisto()
      call vegas1(EvalCS_1L_Zprime_ttbqqb,VG_Result,VG_Error,VG_Chi2)
   endif
ELSEIF( CORRECTION.EQ.4 .AND. PROCESS.EQ.62 ) THEN
   print *, 'Virtual in decay, qqb->Zprime->ttb'
   call vegas(EvalCS_NLODK_Zprime_ttb,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
   itmx = VegasIt1
   ncall= VegasNc1
   call InitHisto()
   call vegas1(EvalCS_NLODK_Zprime_ttb,VG_Result,VG_Error,VG_Chi2)
  endif
ENDIF
ENDIF

IF( MASTERPROCESS.EQ.63 ) THEN
IF( CORRECTION.EQ.2 .AND. (PROCESS.EQ.63 .OR. PROCESS.EQ.64 .OR. PROCESS.EQ.66 )) THEN
   print *, 'Real correction, qqb->Zprime->ttb'
  call vegas(EvalCS_Real_Zprime_ttbqqb,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
    itmx = VegasIt1
    ncall= VegasNc1
    call InitHisto()
    call vegas1(EvalCS_Real_Zprime_ttbqqb,VG_Result,VG_Error,VG_Chi2)
  endif
ENDIF
ENDIF

IF ( MASTERPROCESS.EQ.2) THEN
IF ( CORRECTION.EQ.1 .AND. PROCESS.EQ.65 ) THEN
   call qlinit()
   print *, "Gluon-Z' interference, virtual"
  call vegas(EvalCS_Virt_Zprime_Interf,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
    itmx = VegasIt1
    ncall= VegasNc1
    call InitHisto()
    call vegas1(EvalCS_Virt_Zprime_Interf,VG_Result,VG_Error,VG_Chi2)
  endif
ENDIF
ENDIF

IF ( MASTERPROCESS.EQ.4) THEN
IF ( CORRECTION.EQ.2 .AND. PROCESS.EQ.67 ) THEN
   call qlinit()
   print *, "Gluon-Z' interference, real"
  call vegas(EvalCS_Real_Zprime_Interf,VG_Result,VG_Error,VG_Chi2)
  if( warmup ) then
    itmx = VegasIt1
    ncall= VegasNc1
    call InitHisto()
    call vegas1(EvalCS_Real_Zprime_Interf,VG_Result,VG_Error,VG_Chi2)
  endif
ENDIF
ENDIF


!!! End Zprime section


return
END SUBROUTINE
!DEC$ ENDIF








!DEC$ IF(_UseCuba .EQ.1)
SUBROUTINE StartVegas(VG_Result,VG_Error)
use ModMisc
use ModCrossSection_TTB
use ModCrossSection_TTBJ
use ModCrossSection_TTBP
use ModCrossSection_TTBETmiss
use ModCrossSection_ZprimeTTB
use ModKinematics
use ModParameters
implicit none
include "vegas_common.f"
real(8) :: VG_Result,VG_Error,VG_Chi2
integer :: verbose,mineval,maxeval,nstart,nincrease,neval,fail,ncomp,nbatch,userdata,seed,gridno
real(8) :: epsrel,epsabs
character :: CubaStateFile*(20)



  ncomp = 1
  mineval   = 1000
  maxeval   = VegasNc1
  nstart    = 1000
  nincrease = 1000
  epsrel = 1d-10
  epsabs = 1d-10
  verbose = 2
  userdata = 0
  seed = 0
  nbatch = 10000
  gridno = 1
  cubastatefile = ""
! in bash shell: export CUBACORES=8


  VegasMxDim=mxdim
  call InitHisto()


IF( MASTERPROCESS.EQ.1 ) THEN
IF( CORRECTION   .EQ.0 ) THEN

  call vegas(ndim,ncomp,EvalCS_1L_ttbgg_CUBA,userdata,epsrel,epsabs,verbose,seed,mineval,maxeval,nstart,nincrease,nbatch,gridno,cubastatefile,neval,fail,VG_Result,VG_Error,VG_Chi2)

ENDIF
ENDIF



IF( MASTERPROCESS.EQ.5 ) THEN
IF( CORRECTION   .EQ.2 ) THEN

  call vegas(ndim,ncomp,EvalCS_Real_ttbgggg_CUBA,userdata,epsrel,epsabs,verbose,seed,mineval,maxeval,nstart,nincrease,nbatch,gridno,cubastatefile,neval,fail,VG_Result,VG_Error,VG_Chi2)

ENDIF
ENDIF



IF( MASTERPROCESS.EQ.6 ) THEN
IF( CORRECTION   .EQ.2 ) THEN

  call vegas(ndim,ncomp,EvalCS_Real_ttbqqbgg_CUBA,userdata,epsrel,epsabs,verbose,seed,mineval,maxeval,nstart,nincrease,nbatch,gridno,cubastatefile,neval,fail,VG_Result,VG_Error,VG_Chi2)

ENDIF
ENDIF


it=4


RETURN
END SUBROUTINE
!DEC$ ENDIF






SUBROUTINE InitVegas()
use ModKinematics
use ModParameters
implicit none
include "vegas_common.f"

  VegasIt0_default = 3
  VegasIt1_default = 5

  idum = -VegasSeed
  xl(1:mxdim) = 0d0
  xu(1:mxdim) = 1d0
  acc = -1d0
  nprn = 1
  readin=.false.
  writeout=.false.

  if( VegasIt0.eq.-1 ) VegasIt0 = VegasIt0_default
  if( VegasNc0.eq.-1 ) VegasNc0 = VegasNc0_default
  if( VegasIt1.eq.-1 ) VegasIt1 = VegasIt1_default
  if( VegasNc1.eq.-1 ) VegasNc1 = VegasNc1_default

return
END SUBROUTINE





SUBROUTINE InitPDFs()
use ModParameters
implicit none


IF( PDFSET.EQ.1 ) THEN! MRST/MSTW
! no initialization necessary
  IF( NLOPARAM.EQ.2) THEN
      PDFSetString = " MSTW2008 NLO (mstw2008nlo.00.dat)"
  ELSEIF( NLOPARAM.EQ.0 .OR. NLOPARAM.EQ.1 ) THEN
      PDFSetString = " MSTW2008 LO (mstw2008lo.00.dat)"
  ENDIF



ELSEIF( PDFSET  .EQ.2 ) THEN! CTEQ
  IF( NLOPARAM.EQ.2) THEN
      call SetCtq6(1) !  CTEQ6M   Standard MSbar scheme   0.118     326   226    cteq6m.tbl
!       call SetCtq6(200) !  updated CTEQ6.1M Standard MSbar scheme   0.118     326   226    cteq6m.tbl
!      call SetCtq6(400) !  CTEQ6.6M;                        0.118     326   226    ctq66.00.pds

print *, "SWITCHED TO CTEQ6M FOR CHINESE CHECK"

     PDFSetString = " CTEQ6.6M NLO (ctq66.00.pds)"
!      call SetCT10(100)!   Central CT10           0.118      ct10.00.pds
!      PDFSetString = "CTEQ10 NLO (ct10.00.pds)"
!      print *, "check that in cteq2mrst.f ct10 is really called!"; stop
  ELSEIF( NLOPARAM.EQ.0 .OR. NLOPARAM.EQ.1 ) THEN
     call SetCtq6(4)  ! CTEQ6L1  Leading Order           0.130**   215** 165    cteq6l1.tbl
  !     call SetCtq6(3)  ! CTEQ6L  Leading Order            0.118**   326   226    cteq6l.tbl
     PDFSetString = " CTEQ6L1 LO (cteq6l.tbl)"
  ENDIF
ENDIF


return
END SUBROUTINE




SUBROUTINE OpenFiles()
use ModParameters
implicit none
character :: filename*(100)

!    filename = trim(HistoFile)//'.dat'
!    open(unit=14,file=trim(filename),form='formatted',access= 'sequential',status='replace')            ! Histogram file

   filename = trim(HistoFile)//'.status'
   open(unit=15,file=trim(filename),form='formatted',access= 'sequential',status='replace')         ! Vegas status file

!    filename = trim(HistoFile)//'.tmp_histo'
!    open(unit=16,file=trim(filename),form='formatted',access= 'sequential',status='replace')         ! Histo status file

return
END SUBROUTINE



SUBROUTINE CloseFiles()
implicit none

!    close(14)
   close(15)
!    close(16)

return
END SUBROUTINE




SUBROUTINE WriteHisto(TheUnit,curit,VG_Result,VG_Error,RunTime)
use ModKinematics
use ModParameters
implicit none
include "vegas_common.f"
integer :: NBin,Hits,NHisto,SumHits,TheUnit,curit
real(8) :: BinSize,LowVal,BinVal,Value,Error,Integral
real(8),parameter :: ToGeV=1d2, ToPb=1d-3
real(8) :: VG_Result,VG_Error,RunTime
character :: filename*(100)

  filename = trim(HistoFile)//'.dat'
  if(TheUnit.ne.6) open(unit=TheUnit,file=trim(filename),form='formatted',access= 'sequential',status='replace')   ! Histogram file

  call WriteParameters(TheUnit)
  write(TheUnit,"(A,2X,1F9.2,A)") "# run time =",RunTime/60d0,"min"
  write(TheUnit,"(A,2X,1F20.10)") "# EvalCounter  =",dble(EvalCounter)
  write(TheUnit,"(A,2X,1F20.10)") "# PSCutCounter =",dble(PSCutCounter)
  write(TheUnit,"(A,2X,1F20.10)") "# SkipCounter  =",dble(SkipCounter)

  write(TheUnit,"(A,2X,1PE20.10,2X,1PE20.5)") "#TotCS[fb]=",VG_Result,VG_Error
  do NHisto=1,NumHistograms
      write(TheUnit,"(A,I2,A,A)") "# Histogram",NHisto,": ",Histo(NHisto)%Info
      Integral = 0d0
      SumHits = 0
      BinSize = Histo(NHisto)%BinSize * Histo(NHisto)%SetScale
      LowVal  = Histo(NHisto)%LowVal  * Histo(NHisto)%SetScale
      do NBin=1, Histo(NHisto)%NBins
          BinVal = (LowVal+(NBin-1)*BinSize)
          Hits   = Histo(NHisto)%Hits(NBin)
          SumHits = SumHits + Hits
          if( unweighted ) then
              Value  = Histo(NHisto)%Value(NBin)/BinSize
              Integral = Integral + Histo(NHisto)%Value(NBin)
              Error  = 1d0/dsqrt(dble(Hits))
          else
              Value  = Histo(NHisto)%Value(NBin)/BinSize/curit
              Integral = Integral + Histo(NHisto)%Value(NBin)/curit
              Error  = 1d0/(BinSize)/curit * dsqrt( Histo(NHisto)%Value2(NBin) - 1d0/curit/ncall*Histo(NHisto)%Value(NBin)**2 )
          endif
          if(Hits.ge.999999999) Hits=999999999
          write(TheUnit,"(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)") NHisto,"|",BinVal,"|",Value,"|",Error,"|",Hits,"|"
      enddo
      Integral = Integral + (Histo(NHisto)%Value(0)+Histo(NHisto)%Value(Histo(NHisto)%NBins+1))/curit
      write(TheUnit,"(A,2X,1PE23.16)") "# integrated result:",Integral
      write(TheUnit,"(A,2X,1I23)") "# total number of hits:",SumHits
  enddo

  if(TheUnit.ne.6) close(TheUnit)

return
END SUBROUTINE



SUBROUTINE WriteInstantHisto()
use ModKinematics
use ModParameters
implicit none
include "vegas_common.f"
integer :: NBin,Hits,NHisto
real(8) :: BinSize,LowVal,BinVal,Value,Error,Integral
real(8),parameter :: ToGeV=1d2, ToPb=1d-3
real(8) :: VG_Result,VG_Error,RunTime
character :: filename*(100)

  filename = trim(HistoFile)//'.tmp_histo'
  open(unit=16,file=trim(filename),form='formatted',access= 'sequential',status='replace')         ! Histo status file

  write(16,"(A)") "#"
  write(16,"(A,I3)") "# temporary histogram ",it-1
  write(16,"(A,2X,1F20.10)") "# EvalCounter  =",dble(EvalCounter)
  write(16,"(A,2X,1F20.10)") "# PSCutCounter =",dble(PSCutCounter)
  write(16,"(A,2X,1F20.10)") "# SkipCounter  =",dble(SkipCounter)
  do NHisto=1,NumHistograms
      write(16,"(A,I2,A,A)") "# Histogram",NHisto,": ",Histo(NHisto)%Info
      Integral = 0d0
      BinSize = Histo(NHisto)%BinSize * Histo(NHisto)%SetScale
      LowVal  = Histo(NHisto)%LowVal  * Histo(NHisto)%SetScale
      do NBin=1, Histo(NHisto)%NBins
          BinVal = (LowVal+(NBin-1)*BinSize)
          Hits   = Histo(NHisto)%Hits(NBin)
          Value  = Histo(NHisto)%Value(NBin)/BinSize/(it-1)
          Integral = Integral + Histo(NHisto)%Value(NBin)/(it-1)
          Error  = 1d0/(BinSize)/(it-1) * dsqrt( Histo(NHisto)%Value2(NBin) - 1d0/(it-1)/ncall*Histo(NHisto)%Value(NBin)**2 )
          write(16,"(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)") NHisto,"|",BinVal,"|",Value,"|",Error,"|",Hits,"|"
      enddo
      Integral = Integral + (Histo(NHisto)%Value(0)+Histo(NHisto)%Value(Histo(NHisto)%NBins+1))/(it-1)
      write(16,"(A,2X,1PE23.16)") "# integrated result:",Integral
      write(16,"(A)") "#"
  enddo

  close(16)

return
END SUBROUTINE




SUBROUTINE PlotVegas()
use ModParameters
use ifport
implicit none
include "vegas_common.f"
character :: filename*(100),istr*(2)
integer :: res

    filename = trim(HistoFile)
    write(istr,"(I2)") it

    res = system("./misc/PlotVegasRun.sh "//trim(filename)//" "//trim(istr))
    if( res.eq.-1 ) then
        print *, "Error plotting vegas convergence"
    else
        print *, "Vegas convergence has been plotted in ",trim(filename)//".eps"
    endif


return
END SUBROUTINE


