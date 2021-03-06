# configure generated Makefile
#
# Makefile for Packmol: Read the comments if you have some
#                       problem while compiling.
#
# You may use the ./configure script to search automatically for
# some fortran compiler.
#
# This make file will try to compile packmol with the default
# fortran compiler, defined by the FC directive. For doing this,
# just type
#
#          make 
#
# The default compilation compiles only the serial version of packmol.
#
# If you want to compile with some specific fortran compiler, or
# to enable the parallel implementation, you must change the line
# below to the path of your fortran compiler. The parallel version
# must be compiled with gfortran version 4.2 or greater. 
#
FORTRAN = /usr/bin/gfortran
#
FORTRAN_PARALLEL = gfortran 
#
#
# Change "AUTO" to the fortran command you want. After changing
# this line, you have two options: compile the parallel version
# of packmol (if the compiler is gfortran >= 4.2 or other openmp
# compatible compiler), or compile the serial version.
# To compile the parallel version type
#
#          make parallel
#
#
# Change the flags of the compilation if you want:
#
SERIALFLAGS = -O3 -ffast-math 
#
# Flags for the compilation of the parallel version
#
PARALLELFLAGS = -O3 -ffast-math -fopenmp 
 
###################################################################
#                                                                 #
# Generally no modifications are required after this.             #
#                                                                 #
###################################################################
#
# Get the default fortran compiler
#
ifeq ($(FORTRAN),AUTO)
FORTRAN = $(FC)
endif 
#
# Files required
#
srcall = cenmass.f gencan.f initial.f io.f fgcommon.f packmolnew.f \
         polartocart.f heuristics.f  
srcserial = feasy.f geasy.f
srcparallel = feasyparallel.f geasyparallel.f compindexes.f
commons = Makefile sizes.i gencan.i molpa.i
#
# Compiling the serial version
#
serial : $(srcall) $(srcserial) $(commons)
	@echo " ------------------------------------------------------ " 
	@echo " Compiling packmol with $(FORTRAN) " 
	@echo " Flags: $(SERIALFLAGS) " 
	@echo " ------------------------------------------------------ " 
	@$(FORTRAN) $(SERIALFLAGS) $(srcall) $(srcserial) -o packmol
	@chmod +x solvate.tcl 
	@echo " Packmol succesfully built." 
	@echo " ------------------------------------------------------ " 
#
# Compiling the parallel version: under development.
#
parallel : $(srcall) $(srcparallel) $(commons) ppackmol
	@echo " ------------------------------------------------------ " 
	@echo " Compiling packmol with $(FORTRAN_PARALLEL) " 
	@echo " Flags: $(PARALLELFLAGS)" 
	@echo " ------------------------------------------------------ "
	@$(FORTRAN_PARALLEL) $(PARALLELFLAGS) $(srcall) $(srcparallel)\
                           -o packmol
	@chmod +x ppackmol
	@chmod +x solvate.tcl 
	@echo " ------------------------------------------------------ " 
	@echo " Packmol succesfully built. Paralell version available. " 
	@echo " ------------------------------------------------------ " 
#
# Compiling with flags for development
#
devel : $(srcall) $(srcserial) $(srcparallel) $(commons) ppackmol
	@echo " ------------------------------------------------------ " 
	@echo " Compiling packmol with $(FORTRAN) " 
	@echo " Flags: -Wunused"
	@echo " ------------------------------------------------------ "
	@$(FORTRAN) -Wunused $(srcall) $(srcserial) -o packmol
	@chmod +x solvate.tcl 
	@echo " ------------------------------------------------------ " 
	@echo " Packmol succesfully built. " 
	@echo " ------------------------------------------------------ " 
