# ! TODO: at the end add dependencies for all .o targets according to use statements in modules
Here = $(PWD)
ConfigFile = $(Here)/ttbarjets.cfg
ModuleDir = $(Here)/modules
ObjectDir = $(Here)/objects
DipoleDir = $(Here)/dipoles
OD = $(ObjectDir)
PDFDir = $(Here)/PDFS
VegasDir = $(Here)/Vegas
CubaDir = $(Here)/Cuba
OptReport = $(Here)/OptRep.txt
PSDir = $(Here)/PhaseSpace
QCDLoop = $(Here)/QCDLoop-1.9


# compiler options
F95compiler = ifort
Opt = Yes


ifeq ($(Opt),Yes)
   IfortOpts   = -O1 -fpp -vec-report0 -opt-report -opt-report-file$(OptReport) -I$(Here)/colors -I$(VegasDir) -module $(ModuleDir)
else
   IfortOpts   = -O0 -fpp -implicitnone -zero -check bounds -check pointer -warn interfaces -ftrapuv -I$(Here)/colors -I$(VegasDir) -module $(ModuleDir)
endif
fcomp = $(F95compiler) $(IfortOpts) @$(ConfigFile)

# never use gcc with other than O0, optimization seems to be buggy
ccomp = gcc -O0


# executable file
Exec = ./TOPAZ


makeDep = $(ConfigFile) \
          makefile


# fastjet stuff
#FASTJET_CONFIG=/home/schulze/usr/local/bin/fastjet-config

FASTJET_CONFIG=/home/raoul/Work/fastjet-install/bin/fastjet-config
CXXFLAGS += $(shell $(FASTJET_CONFIG) --cxxflags)
FJLIBS += $(shell $(FASTJET_CONFIG) --libs --plugins )


PDFObj = $(PDFDir)/mstwpdf.o \
         $(PDFDir)/mrst2001lo.o \
         $(PDFDir)/cteq2mrst.o \
         $(PDFDir)/Cteq66Pdf.o \
         $(PDFDir)/CT10Pdf.o
RockyObj = $(PSDir)/genps.o \
           $(PSDir)/boost.o
YetiObj  = $(PSDir)/yeti.o
VegasObj = $(VegasDir)/vegas.o

CubaLib = $(CubaDir)/libcuba.a



IntegralObj = $(QCDLoop)/ql/libqcdloop.a\
              $(QCDLoop)/ff/libff.a


OPPDep = mod_NVBasis.f90 \
         mod_Residues.f90 \
         mod_UCuts.f90 \
         mod_Residues_new.f90 \
         mod_UCuts_new.f90

OPPObj = $(ObjectDir)/mod_NVBasis.o \
         $(ObjectDir)/mod_Residues.o \
         $(ObjectDir)/mod_UCuts.o \
         $(ObjectDir)/mod_Residues_new.o \
         $(ObjectDir)/mod_UCuts_new.o \
         $(ObjectDir)/mod_NVBasis128.o \
         $(ObjectDir)/mod_Residues128.o \
         $(ObjectDir)/mod_UCuts128.o \
         $(ObjectDir)/mod_Residues128_new.o \
         $(ObjectDir)/mod_UCuts128_new.o 

DipoleDepTTB = $(Here)/dipoles/mod_Dipoles_GGTTBG.f90 \
            $(Here)/dipoles/mod_Dipoles_GGTTBG_noDK.f90 \
            $(Here)/dipoles/mod_Dipoles_QGTTBQ.f90 \
            $(Here)/dipoles/mod_Dipoles_QGTTBQ_noDK.f90 \
            $(Here)/dipoles/mod_Dipoles_QQBTTBG.f90 \
            $(Here)/dipoles/mod_Dipoles_QQBTTBG_noDK.f90 \
            $(Here)/dipoles/mod_IntDipoles_GGTTBG.f90 \
            $(Here)/dipoles/mod_IntDipoles_GGTTBG_noDK.f90 \
            $(Here)/dipoles/mod_IntDipoles_QGTTBQ.f90 \
            $(Here)/dipoles/mod_IntDipoles_QGTTBQ_noDK.f90 \
            $(Here)/dipoles/mod_IntDipoles_QQBTTBG.f90 \
            $(Here)/dipoles/mod_IntDipoles_QQBTTBG_noDK.f90

DipoleObjTTB = $(ObjectDir)/mod_Dipoles_GGTTBG.o \
            $(ObjectDir)/mod_Dipoles_GGTTBG_noDK.o \
            $(ObjectDir)/mod_Dipoles_QGTTBQ.o \
            $(ObjectDir)/mod_Dipoles_QGTTBQ_noDK.o \
            $(ObjectDir)/mod_Dipoles_QQBTTBG.o \
            $(ObjectDir)/mod_Dipoles_QQBTTBG_noDK.o \
            $(ObjectDir)/mod_IntDipoles_GGTTBG.o \
            $(ObjectDir)/mod_IntDipoles_GGTTBG_noDK.o \
            $(ObjectDir)/mod_IntDipoles_QGTTBQ.o \
            $(ObjectDir)/mod_IntDipoles_QGTTBQ_noDK.o \
            $(ObjectDir)/mod_IntDipoles_QQBTTBG.o \
            $(ObjectDir)/mod_IntDipoles_QQBTTBG_noDK.o




DipoleDepTTBJ = $(Here)/dipoles/mod_Dipoles_GGTTBGG2.f90 \
            $(Here)/dipoles/mod_Dipoles_GGTTBQQB2.f90 \
            $(Here)/dipoles/mod_Dipoles_QQBTTBGG2.f90 \
            $(Here)/dipoles/mod_Dipoles_QGTTBQG2.f90 \
            $(Here)/dipoles/mod_Dipoles_QBGTTBQBG2.f90 \
            $(Here)/dipoles/mod_Dipoles_qbb_ttqqqq_dip2.f90 \
            $(Here)/dipoles/mod_Dipoles_qqb_ttqqqq_dip2.f90 \
            $(Here)/dipoles/mod_Dipoles_qqq_ttqqqq_dip2.f90 \
            $(Here)/dipoles/mod_Dipoles_DKJ_TTB.f90 \
            $(Here)/dipoles/mod_Dipoles_DKJ_GGTTBG.f90 \
            $(Here)/dipoles/mod_Dipoles_DKJ_QQBTTBG.f90 \
            $(Here)/dipoles/mod_Dipoles_DKJ_QGTTBQ.f90 \
            $(Here)/dipoles/mod_IntDipoles_GGTTBGG2.f90 \
            $(Here)/dipoles/mod_IntDipoles_GGTTBQQB2.f90 \
            $(Here)/dipoles/mod_IntDipoles_QQBTTBGG2.f90 \
            $(Here)/dipoles/mod_IntDipoles_QGTTBQG2.f90 \
            $(Here)/dipoles/mod_IntDipoles_QBGTTBQBG2.f90 \
            $(Here)/dipoles/mod_IntDipoles_SixQuark2.f90 \
            $(Here)/dipoles/mod_IntDipoles_DKJ_GGTTBG.f90 \
            $(Here)/dipoles/mod_IntDipoles_DKJ_QQBTTBG.f90 \
            $(Here)/dipoles/mod_IntDipoles_DKJ_QGTTBQ.f90


DipoleObjTTBJ = $(ObjectDir)/mod_Dipoles_GGTTBGG2.o \
            $(ObjectDir)/mod_Dipoles_GGTTBQQB2.o \
            $(ObjectDir)/mod_Dipoles_QQBTTBGG2.o \
            $(ObjectDir)/mod_Dipoles_QGTTBQG2.o \
            $(ObjectDir)/mod_Dipoles_QBGTTBQBG2.o \
            $(ObjectDir)/mod_Dipoles_qbb_ttqqqq_dip2.o \
            $(ObjectDir)/mod_Dipoles_qqb_ttqqqq_dip2.o \
            $(ObjectDir)/mod_Dipoles_qqq_ttqqqq_dip2.o \
            $(ObjectDir)/mod_Dipoles_DKJ_TTB.o \
            $(ObjectDir)/mod_Dipoles_DKJ_GGTTBG.o \
            $(ObjectDir)/mod_Dipoles_DKJ_QQBTTBG.o \
            $(ObjectDir)/mod_Dipoles_DKJ_QGTTBQ.o \
            $(ObjectDir)/mod_IntDipoles_GGTTBGG2.o \
            $(ObjectDir)/mod_IntDipoles_GGTTBQQB2.o \
            $(ObjectDir)/mod_IntDipoles_QQBTTBGG2.o \
            $(ObjectDir)/mod_IntDipoles_QGTTBQG2.o \
            $(ObjectDir)/mod_IntDipoles_QBGTTBQBG2.o \
            $(ObjectDir)/mod_IntDipoles_SixQuark2.o \
            $(ObjectDir)/mod_IntDipoles_DKJ_GGTTBG.o \
            $(ObjectDir)/mod_IntDipoles_DKJ_QQBTTBG.o \
            $(ObjectDir)/mod_IntDipoles_DKJ_QGTTBQ.o





DipoleDepTTBP = $(Here)/dipoles/mod_Dipoles_GGTTBGP.f90 \
            $(Here)/dipoles/mod_Dipoles_QQBTTBGP.f90 \
            $(Here)/dipoles/mod_Dipoles_QGTTBQP.f90 \
            $(Here)/dipoles/mod_Dipoles_QBGTTBQBP.f90 \
            $(Here)/dipoles/mod_Dipoles_DKP_GGTTBG.f90 \
            $(Here)/dipoles/mod_Dipoles_DKP_QQBTTBG.f90 \
            $(Here)/dipoles/mod_Dipoles_DKP_QGTTBQ.f90 \
            $(Here)/dipoles/mod_IntDipoles_GGTTBGP.f90 \
            $(Here)/dipoles/mod_IntDipoles_QQBTTBGP.f90 \
            $(Here)/dipoles/mod_IntDipoles_QGTTBQP.f90 \
            $(Here)/dipoles/mod_IntDipoles_QBGTTBQBP.f90 \
            $(Here)/dipoles/mod_IntDipoles_DKP_GGTTBG.f90 \
            $(Here)/dipoles/mod_IntDipoles_DKP_QQBTTBG.f90 \
            $(Here)/dipoles/mod_IntDipoles_DKP_QGTTBQ.f90

DipoleObjTTBP = $(ObjectDir)/mod_Dipoles_GGTTBGP.o \
                $(ObjectDir)/mod_Dipoles_QQBTTBGP.o \
                $(ObjectDir)/mod_Dipoles_QGTTBQP.o \
                $(ObjectDir)/mod_Dipoles_QBGTTBQBP.o \
                $(ObjectDir)/mod_Dipoles_DKP_GGTTBG.o \
                $(ObjectDir)/mod_Dipoles_DKP_QQBTTBG.o \
                $(ObjectDir)/mod_Dipoles_DKP_QGTTBQ.o \
                $(ObjectDir)/mod_IntDipoles_GGTTBGP.o \
                $(ObjectDir)/mod_IntDipoles_QQBTTBGP.o \
                $(ObjectDir)/mod_IntDipoles_QGTTBQP.o \
                $(ObjectDir)/mod_IntDipoles_QBGTTBQBP.o \
                $(ObjectDir)/mod_IntDipoles_DKP_GGTTBG.o \
                $(ObjectDir)/mod_IntDipoles_DKP_QQBTTBG.o \
                $(ObjectDir)/mod_IntDipoles_DKP_QGTTBQ.o



DipoleDepTTBZ = $(Here)/dipoles/mod_Dipoles_GGTTBGZ.f90 \
            $(Here)/dipoles/mod_Dipoles_QQBTTBGZ.f90 \
            $(Here)/dipoles/mod_Dipoles_QGTTBQZ.f90 \
            $(Here)/dipoles/mod_Dipoles_QBGTTBQBZ.f90 \
            $(Here)/dipoles/mod_IntDipoles_GGTTBGZ.f90 \
            $(Here)/dipoles/mod_IntDipoles_QQBTTBGZ.f90 \
            $(Here)/dipoles/mod_IntDipoles_QGTTBQZ.f90 \
            $(Here)/dipoles/mod_IntDipoles_QBGTTBQBZ.f90

DipoleObjTTBZ = $(ObjectDir)/mod_Dipoles_GGTTBGZ.o \
                $(ObjectDir)/mod_Dipoles_QQBTTBGZ.o \
                $(ObjectDir)/mod_Dipoles_QGTTBQZ.o \
                $(ObjectDir)/mod_Dipoles_QBGTTBQBZ.o \
                $(ObjectDir)/mod_IntDipoles_GGTTBGZ.o \
                $(ObjectDir)/mod_IntDipoles_QQBTTBGZ.o \
                $(ObjectDir)/mod_IntDipoles_QGTTBQZ.o \
                $(ObjectDir)/mod_IntDipoles_QBGTTBQBZ.o


DipoleDepSTSTB = $(Here)/dipoles/mod_Dipoles_GGSTSTBG.f90 \
		 $(Here)/dipoles/mod_IntDipoles_GGSTSTBG.f90 \
		 $(Here)/dipoles/mod_Dipoles_QQBSTSTBG.f90  \
		 $(Here)/dipoles/mod_Dipoles_QGSTSTBQ.f90 \
                 $(Here)/dipoles/mod_IntDipoles_QQBSTSTBG.f90 \
                 $(Here)/dipoles/mod_IntDipoles_QGSTSTBQ.f90

DipoleObjSTSTB = $(ObjectDir)/mod_Dipoles_GGSTSTBG.o \
		 $(ObjectDir)/mod_IntDipoles_GGSTSTBG.o \
		 $(ObjectDir)/mod_Dipoles_QQBSTSTBG.o  \
		 $(ObjectDir)/mod_Dipoles_QGSTSTBQ.o \
                 $(ObjectDir)/mod_IntDipoles_QQBSTSTBG.o \
                 $(ObjectDir)/mod_IntDipoles_QGSTSTBQ.o




DipoleDepHTHTB = $(Here)/dipoles/mod_Dipoles_GGHTHTBG.f90 \
		 $(Here)/dipoles/mod_Dipoles_QQBHTHTBG.f90  \
		 $(Here)/dipoles/mod_Dipoles_QGHTHTBQ.f90  \
		 $(Here)/dipoles/mod_IntDipoles_GGHTHTBG.f90 \
                 $(Here)/dipoles/mod_IntDipoles_QQBHTHTBG.f90 \
                 $(Here)/dipoles/mod_IntDipoles_QGHTHTBQ.f90

DipoleObjHTHTB = $(ObjectDir)/mod_Dipoles_GGHTHTBG.o \
		 $(ObjectDir)/mod_Dipoles_QQBHTHTBG.o  \
		 $(ObjectDir)/mod_Dipoles_QGHTHTBQ.o \
		 $(ObjectDir)/mod_IntDipoles_GGHTHTBG.o \
                 $(ObjectDir)/mod_IntDipoles_QQBHTHTBG.o \
                 $(ObjectDir)/mod_IntDipoles_QGHTHTBQ.o



# Zprime section

DipoleDepZprime = $(Here)/dipoles/mod_Dipoles_ZprimeTTB.f90 \
		  $(Here)/dipoles/mod_IntDipoles_ZprimeTTB.f90

DipoleObjZprime = $(ObjectDir)/mod_Dipoles_ZprimeTTB.o \
		 $(ObjectDir)/mod_IntDipoles_ZprimeTTB.o 

# End Zprime section


MadGraphObj = $(Here)/MadGraph/gg_ttbg.o \
				  $(Here)/MadGraph/ddb_Zprime_ttb.o \
				  $(Here)/MadGraph/uub_Zprime_ttb.o \
				  $(Here)/MadGraph/switchmom.o \
				  $(HOME)/Lib/HELAS-3.0/coupsm.o \
				  $(HOME)/Lib/HELAS-3.0/oxxxxx.o \
				  $(HOME)/Lib/HELAS-3.0/vxxxxx.o \
				  $(HOME)/Lib/HELAS-3.0/ixxxxx.o \
				  $(HOME)/Lib/HELAS-3.0/fvixxx.o \
				  $(HOME)/Lib/HELAS-3.0/jvvxxx.o \
				  $(HOME)/Lib/HELAS-3.0/jgggxx.o \
				  $(HOME)/Lib/HELAS-3.0/jioxxx.o \
				  $(HOME)/Lib/HELAS-3.0/iovxxx.o \
				  $(HOME)/Lib/HELAS-3.0/fvoxxx.o \
				  $(HOME)/Lib/HELAS-3.0/vvvxxx.o \
				  $(HOME)/Lib/HELAS-3.0/libdhelas3.ifc90.a \
				  $(Here)/MadGraph/gg_ttbz.o \
				  $(Here)/MadGraph/uub_ttbz.o \
				  $(Here)/MadGraph/ddb_ttbz.o
#				  $(Here)/MadGraph/gg_ttb.o \
#				  $(Here)/MadGraph/gg_ttba.o \
#				  $(Here)/MadGraph/gg_tbtga.o \
#				  $(Here)/MadGraph/ug_ttbua.o \
#				  $(Here)/MadGraph/ddb_ttba.o \
#				  $(Here)/MadGraph/uub_ttba.o \
#				  $(Here)/MadGraph/dbg_ttbdba.o \
#				  $(Here)/MadGraph/ddb_ttbga.o \
#				  $(Here)/MadGraph/wm_ubduub.o  \
#				  $(Here)/MadGraph/wm_ubdccb.o  \
#				  $(Here)/MadGraph/gg_ttbgg.o \
#				  $(Here)/MadGraph/tb_bbemveb.o \
#				  $(Here)/MadGraph/tb_tbepem.o \
#				  $(Here)/MadGraph/gg_bepvebbemve.o \
#				  $(Here)/MadGraph/gg_uub.o \
#				  $(Here)/MadGraph/gg_uubg.o \
#				  $(Here)/MadGraph/uub_ttb.o \
#				  $(Here)/MadGraph/uub_ttbg.o \
#				  $(Here)/MadGraph/uub_ttbgg.o \
#				  $(Here)/MadGraph/uub_ttbddb.o \
#				  $(Here)/MadGraph/uub_ttbuub.o \
#				  $(Here)/MadGraph/ubdb_ttbubdb.o \
#				  $(Here)/MadGraph/udb_ttbudb.o \
#				  $(Here)/MadGraph/ud_ttbud.o \
#				  $(Here)/MadGraph/uu_ttbuu.o \
#				  $(Here)/MadGraph/ubub_ttbubub.o \
#				  $(Here)/MadGraph/gg_ttbuub.o \
#				  $(Here)/MadGraph/uub_ddbg.o \
#				  $(Here)/MadGraph/ug_ttbu.o \
#				  $(Here)/MadGraph/ubg_ttbub.o \
#				  $(Here)/MadGraph/gub_ttbub.o \
#				  $(Here)/MadGraph/ug_ttbug.o \
#				  $(Here)/MadGraph/gu_ttbu.o \
#				  $(Here)/MadGraph/ubg_ttbubg.o \
#				  $(Here)/MadGraph/t_bepve.o \
#				  $(Here)/MadGraph/tb_emvebbb.o \
#				  $(Here)/MadGraph/tb_bbemveb.o \
#				  $(Here)/MadGraph/tb_bbemvebgg.o \
#				  $(Here)/MadGraph/t_bepvegg.o \
#				  $(Here)/MadGraph/t_epvebbbb.o \
#				  $(Here)/MadGraph/t_epvebuub.o \
#				  $(Here)/MadGraph/tb_emvebbbuub.o \
#				  $(Here)/MadGraph/tb_emvebbbbbb.o \
#				  $(Here)/MadGraph/t_bdbuccb.o \
#				  $(Here)/MadGraph/t_bdbuddb.o \
#				  $(Here)/MadGraph/t_bdbuuub.o \
#				  $(Here)/MadGraph/tb_bbubdccb.o \
#				  $(Here)/MadGraph/tb_bbubdddb.o \
#				  $(Here)/MadGraph/tb_bbubduub.o \
#				  $(Here)/MadGraph/tb_bbubdgg.o \
#				  $(Here)/MadGraph/t_bdbugg.o \
#				  $(Here)/MadGraph/tb_bbwm.o  \
#				  $(Here)/MadGraph/t_bwp.o \
#				  $(Here)/MadGraph/t_epveb.o \
#				  $(Here)/MadGraph/t_bwpuub.o \
#				  $(Here)/MadGraph/t_bwpbbb.o \
#				  $(Here)/MadGraph/wm_emveb.o \
#				  $(Here)/MadGraph/wp_epve.o \


# ------------------------------------------------------------


# the order of these object files corresponds to their mutual dependency
allObjects =   				$(ObjectDir)/mod_Misc.o \
					$(ObjectDir)/mod_Parameters.o \
					$(ObjectDir)/mod_Process.o \
					$(ObjectDir)/mod_Permutations.o \
					$(ObjectDir)/mod_IntegerPartition.o \
					$(ObjectDir)/mod_MyRecurrence.o \
					$(ObjectDir)/mod_MyWeylRecurrence.o \
					$(ObjectDir)/mod_Amplitudes.o \
					$(ObjectDir)/mod_Amplitudes_Zprime.o \
					$(OPPObj) \
					$(ObjectDir)/mod_DKIntDipoles.o \
					$(ObjectDir)/mod_HadrWDecay.o \
					$(ObjectDir)/mod_JPsiFrag.o \
					$(Here)/includes_DKP/DKP1L.o \
					$(Here)/includes_DKJ/TopCoeffsNLO.o \
					$(Here)/includes_DKJ/ATopCoeffsNLO.o \
					$(ObjectDir)/mod_TTBP_NLODK.o \
					$(ObjectDir)/mod_TTBJ_NLODK.o \
					$(ObjectDir)/mod_TTBJ_NLODKW.o \
					$(ObjectDir)/mod_WDecay.o \
					$(ObjectDir)/mod_ZDecay.o \
					$(ObjectDir)/mod_TopDecay.o \
					$(ObjectDir)/mod_ExoticDecay.o \
					$(ObjectDir)/mod_IntDipoles.o \
					$(ObjectDir)/mod_SpinCorrel.o \
					$(ObjectDir)/mod_Kinematics.o \
					$(ObjectDir)/mod_Integrals.o \
					$(DipoleObjTTB) \
					$(DipoleObjTTBJ) \
					$(DipoleObjTTBP) \
					$(DipoleObjTTBZ) \
					$(DipoleObjSTSTB) \
					$(DipoleObjHTHTB) \
					$(DipoleObjZprime) \
					$(ObjectDir)/mod_SixFermionProcs2.o \
					$(ObjectDir)/mod_CrossSection_TTB.o \
					$(ObjectDir)/mod_CrossSection_TTBJ.o \
					$(ObjectDir)/mod_CrossSection_TTBP.o \
					$(ObjectDir)/mod_CrossSection_TTBETmiss.o \
					$(ObjectDir)/mod_CrossSection_ZprimeTTB.o \
					$(ObjectDir)/mod_CrossSection_TTBZ.o \
					$(ObjectDir)/main.o


#--------------------------------------------------------------------------------------------------
# note that DKP1L.o,TopCoeffsNLO.o,ATopCoeffsNLO.o are stored in the "include..." directories
# they are not being removed by "make clean"
# their compilation takes long so it should be done only once
#
# external libraries for PDFs, PS Generators and Integrals are not compiled with this makefile
#--------------------------------------------------------------------------------------------------



all: $(allObjects)
	@echo " linking"
	@echo " executable file is " $(Exec)
	@echo " "
# 	$(fcomp) -o $(Exec) $(allObjects) $(RockyObj) $(YetiObj) $(IntegralObj) $(CubaLib) $(PDFObj) 
	$(fcomp) -o $(Exec) $(allObjects) $(RockyObj) $(YetiObj) $(IntegralObj) $(VegasObj) $(PDFObj) #$(MadGraphObj)
# $(ObjectDir)/fastjetfortran.o $(FJLIBS) -lstdc++    add this to above line when fastjet routines are used


clean:
	rm -f ./modules/*.mod
	rm -f ./objects/*.o


# ------------------------------------------------------------





$(ObjectDir)/mod_Misc.o: mod_Misc.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_Parameters.o: mod_Parameters.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_Process.o: mod_Process.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_Permutations.o: mod_Permutations.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_IntegerPartition.o: mod_IntegerPartition.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_MyRecurrence.o: mod_MyRecurrence.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_MyWeylRecurrence.o: mod_MyWeylRecurrence.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_Amplitudes.o: mod_Amplitudes.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@

$(ObjectDir)/mod_Amplitudes_Zprime.o: mod_Amplitudes_Zprime.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@

$(ObjectDir)/main.o: main.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@



$(ObjectDir)/mod_CrossSection_TTB.o: mod_CrossSection_TTB.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_CrossSection_TTBJ.o: mod_CrossSection_TTBJ.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_CrossSection_TTBP.o: mod_CrossSection_TTBP.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_CrossSection_TTBETmiss.o: mod_CrossSection_TTBETmiss.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_CrossSection_ZprimeTTB.o: mod_CrossSection_ZprimeTTB.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_CrossSection_TTBZ.o: mod_CrossSection_TTBZ.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_DKIntDipoles.o: mod_DKIntDipoles.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@



$(ObjectDir)/mod_HadrWDecay.o: mod_HadrWDecay.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_IntDipoles.o: mod_IntDipoles.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@



$(ObjectDir)/mod_Integrals.o: mod_Integrals.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_JPsiFrag.o: mod_JPsiFrag.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_Kinematics.o: mod_Kinematics.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@




$(ObjectDir)/mod_SixFermionProcs2.o: mod_SixFermionProcs2.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_ExoticDecay.o: mod_ExoticDecay.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@



$(ObjectDir)/mod_SpinCorrel.o: mod_SpinCorrel.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_TopDecay.o: mod_TopDecay.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(Here)/includes_DKP/DKP1L.o: $(Here)/includes_DKP/DKP1L.f
	@echo " compiling" $<
	@echo " this typically takes 10 minutes"
	$(fcomp) -80 -I$(Here)/includes_DKP -module $(ModuleDir) -c $< -o $@


$(Here)/includes_DKJ/TopCoeffsNLO.o: $(Here)/includes_DKJ/TopCoeffsNLO.f90
	@echo " compiling" $<
	@echo " this typically takes 10 minutes"
	$(fcomp) -I$(Here)/includes_DKJ -module $(ModuleDir) -c $< -o $@


$(Here)/includes_DKJ/ATopCoeffsNLO.o: $(Here)/includes_DKJ/ATopCoeffsNLO.f90
	@echo " compiling" $<
	@echo " this typically takes 10 minutes"
	$(fcomp) -I$(Here)/includes_DKJ -module $(ModuleDir) -c $< -o $@


$(ObjectDir)/mod_TTBJ_NLODK.o: mod_TTBJ_NLODK.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_TTBJ_NLODKW.o: mod_TTBJ_NLODKW.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@



$(ObjectDir)/mod_TTBP_NLODK.o: mod_TTBP_NLODK.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c -fixed -80 -I$(Here)/includes_DKP $< -o $@



$(ObjectDir)/mod_WDecay.o: mod_WDecay.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@



$(ObjectDir)/mod_ZDecay.o: mod_ZDecay.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@



$(OPPObj): $(OPPDep) $(makeDep)
	fpp @quadprec.cfg mod_NVBasis.f90  mod_NVBasis128.f90
	fpp @quadprec.cfg mod_Residues.f90 mod_Residues128.f90
	fpp @quadprec.cfg mod_UCuts.f90    mod_UCuts128.f90
	fpp @quadprec.cfg mod_Residues_new.f90  mod_Residues128_new.f90
	fpp @quadprec.cfg mod_UCuts_new.f90     mod_UCuts128_new.f90
	@echo " compiling OPP (64bit) "
	$(fcomp) -c mod_NVBasis.f90 -o $(ObjectDir)/mod_NVBasis.o
	$(fcomp) -c mod_Residues.f90 -o $(ObjectDir)/mod_Residues.o
	$(fcomp) -c mod_Residues_new.f90 -o $(ObjectDir)/mod_Residues_new.o
	$(fcomp) -c mod_UCuts.f90 -o $(ObjectDir)/mod_UCuts.o
	$(fcomp) -c mod_UCuts_new.f90 -o $(ObjectDir)/mod_UCuts_new.o
	@echo " compiling OPP (128bit) "
	$(fcomp) -r16 -c mod_NVBasis128.f90 -o $(ObjectDir)/mod_NVBasis128.o
	$(fcomp) -r16 -c mod_Residues128.f90 -o $(ObjectDir)/mod_Residues128.o
	$(fcomp) -r16 -c mod_Residues128_new.f90 -o $(ObjectDir)/mod_Residues128_new.o
	$(fcomp) -r16 -c mod_UCuts128.f90 -o $(ObjectDir)/mod_UCuts128.o
	$(fcomp) -r16 -c mod_UCuts128_new.f90 -o $(ObjectDir)/mod_UCuts128_new.o
	rm mod_NVBasis128.f90 mod_Residues128.f90 mod_UCuts128.f90 mod_UCuts128_new.f90 mod_Residues128_new.f90



$(ObjectDir)/mod_Dipoles%.o: $(DipoleDir)/mod_Dipoles%.f90 $(makeDep)
	@echo " compiling dipoles" mod_Dipoles$*".f90"
	$(fcomp) -c $(DipoleDir)/mod_Dipoles$*".f90" -o $(ObjectDir)/mod_Dipoles$*".o"


$(ObjectDir)/mod_IntDipoles%.o: $(DipoleDir)/mod_IntDipoles%.f90 $(makeDep)
	@echo " compiling dipoles" mod_IntDipoles$*".f90"
	$(fcomp) -c $(DipoleDir)/mod_IntDipoles$*".f90" -o $(ObjectDir)/mod_IntDipoles$*".o"








# $(RockyObj): $(RockyDep) $(makeDep)
# 	@echo " compiling" $< $@
# 	cd $(RockyDir)
# 	$(ccomp) -c $(RockyDep)
# 	cd $(Here)
#
#
# $(YetiObj): $(YetiDep) $(makeDep)
# 	@echo " compiling" $< $@
# 	cd $(YetiDir)
# 	$(fcomp) -c $(YetiDep)
# 	cd $(Here)
#
#
# $(VegasObj): $(VegasDep) $(makeDep)
# 	@echo " compiling" $< $@
# 	cd $(VegasDir)
# 	$(fcomp) -c $(VegasDep)
# 	cd $(Here)
#
#
#
# # pdfs
# $(PDFObj): $(PDFDep)
# 	@echo " compiling pdfs"
# 	cd $(PDFDir)
# 	$(fcomp) -c $(PDFDep)
# 	cd $(Here)



# $(ObjectDir)/genps.o: $(RockyDir)/genps.c $(makeDep)
# 	@echo " compiling" $< $@
# 	$(ccomp) -c $< -o $@
#
#
# $(ObjectDir)/boost.o: $(RockyDir)/boost.c $(makeDep)
# 	@echo " compiling" $< $@
# 	$(ccomp) -c $< -o $@
#
#
# $(ObjectDir)/yeti.o: $(YetiDir)/yeti.f $(makeDep)
# 	@echo " compiling" $<
# 	$(fcomp) -c $< -o $@
#
#
# $(ObjectDir)/vegas.o: $(VegasDep)
# 	@echo " compiling" $<
# 	$(fcomp) -c $< -o $@
#
#
# # pdfs
# $(PDFObj): $(PDFDep)
# 	@echo " compiling pdfs"
# 	$(fcomp) -c $(PDFDep)
#
#

# supresses command calls
.SILENT:
