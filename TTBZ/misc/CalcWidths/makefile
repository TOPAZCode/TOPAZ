Here = $(PWD)
VegasDir = /home/schulze/lib/Vegas2
PSGenDir = $(HOME)/lib/Rocky
ModuleDir = $(Here)/modules
ObjectDir = $(Here)/objects
TTBJDir = $(HOME)/1ttbjets
MakeFile = $(Here)/makefile

# compiler options
F95compiler = ifort
Opt   = Yes


ifeq ($(Opt),Yes)
   IfortOpts   = -O2 -fpp -vec-report0 -I$(VegasDir) -module $(ModuleDir)
else
   IfortOpts   = -O0 -fpp -implicitnone -zero -check bounds -check pointer -warn interfaces -ftrapuv -I$(VegasDir) -module $(ModuleDir)
endif
fcomp = $(F95compiler) $(IfortOpts)

# never use gcc with other than O0, optimization is buggy
ccomp = gcc -O0

# executable file
Exec = ./runInt






# main dependencies & objects
MainDep = mod_Misc.f90 \
          mod_Parameters.f90 \
          mod_Kinematics.f90 \
          mod_Integrand.f90 \
          main.f90

MainObj = $(ObjectDir)/mod_Misc.o \
          $(ObjectDir)/mod_Parameters.o \
          $(ObjectDir)/mod_Kinematics.o \
          $(ObjectDir)/mod_Integrand.o \
          $(ObjectDir)/main.o


MadGraphObj = $(TTBJDir)/MadGraph/t_bwp.o \
              $(TTBJDir)/MadGraph/wp_epve.o \
              $(TTBJDir)/MadGraph/t_bwpg.o \
              $(TTBJDir)/MadGraph/wp_udbg.o \
              $(TTBJDir)/MadGraph/wp_udb.o \
              $(TTBJDir)/MadGraph/switchmom.o \
              $(HOME)/lib/HELAS-3.0/libdhelas3.ifc90.a

# vegas dependencies & objects
VegasDep = $(VegasDir)/vegas.f
VegasObj = $(ObjectDir)/vegas.o


# phase space dependencies & objects
RockyDep = $(PSGenDir)/genps.c \
           $(PSGenDir)/boost.c
RockyObj = $(ObjectDir)/genps.o \
           $(ObjectDir)/boost.o






runInt: $(MakeFile) $(RockyObj) $(VegasObj) $(MainObj)
	@echo " linking"
	@echo " executable file is " $(Exec)
	@echo " "
	$(fcomp) -o  $(Exec) $(MainObj) $(RockyObj) $(VegasObj) $(MadGraphObj)





$(MainObj): $(MakeFile) $(MainDep)
	@echo " compiling runInt"
	$(fcomp) -c $(MainDep)
	mv *.o $(ObjectDir)



# phase space generator
$(RockyObj): $(MakeFile) $(RockyDep)
	@echo " compiling phase space generator and boost routine"
	$(ccomp) -c $(RockyDep)
	mv *.o $(ObjectDir)


# vegas
$(VegasObj): $(MakeFile) $(VegasDep)
	@echo " compiling vegas"
	$(fcomp) -c $(VegasDep)
	mv *.o $(ObjectDir)





# supresses command calls
.SILENT:

