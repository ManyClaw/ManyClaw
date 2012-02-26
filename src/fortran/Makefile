#
# Fortran Makefile for ManyClaw
#

SHELL = /bin/sh
EXECUTABLE = manyclaw_fortran

# Compiler and linkers
FC ?= gfortran
LINK ?= $(FC)
LIB_PATHS ?= 

# Flags, remove ? to override environment settings
FFLAGS ?=
LFLAGS ?= $(FFLAGS)
INCLUDE ?=

# Built flag variables
ALL_INCLUDE ?=
ALL_FFLAGS ?=
ALL_LFLAGS ?=

# Add includes
ALL_INCLUDE += $(INCLUDE)
ALL_INCLUDE += $(addprefix -I,$(MODULE_PATHS)) $(addprefix -I,$(LIB_PATHS))

# Module flag setting, please add other compilers here as necessary
ifeq ($(findstring gfortran,$(FC)),gfortran)
	# There should be no space between this flag and the argument
	MODULE_FLAG = -J
	PPFLAGS = -cpp
    OPT_FLAGS = -O3 -funroll-loops -finline-functions -fdefault-real-8
    DEBUG_FLAGS = -W -Wall -fbounds-check -pedantic-errors -Wunderflow -ffpe-trap=zero,overflow -g
	OMP_FLAG = -fopenmp
else ifeq ($(FC),ifort)
	# Note that there shoud be a space after this flag
	MODULE_FLAG = -module 
	PPFLAGS = -fpp
    OPT_FLAGS = -O2 -ip
    DEBUG_FLAGS = -g -check all -traceback -fpe0
	OMP_FLAGS = -openmp
endif
MAKELEVEL ?= 0
ifeq ($(MAKELEVEL),0)
	ALL_FFLAGS += $(PPFLAGS)
	# ALL_LFLAGS +=
endif
	
# ALL_FFLAGS and ALL_LFLAGS currently only includes the user defined flags
ALL_FFLAGS += $(FFLAGS)
ALL_LFLAGS += $(LFLAGS)

# Source and object lists
MODULES = \
	./precision_module.f90 \
	./test_functions_module.f90 \
	./timer_module.f90

SOURCES = \
	./main.f90

OBJECTS = $(subst .F,.o, $(subst .F90,.o, $(subst .f,.o, $(subst .f90,.o, $(SOURCES)))))
MODULE_FILES = $(subst .F,.mod, $(subst .F90,.mod, $(subst .f,.mod, $(subst .f90,.mod, $(MODULES)))))
MODULE_PATHS = $(sort $(dir $(MODULE_FILES)))
MODULE_OBJECTS = $(subst .F,.o, $(subst .F90,.o, $(subst .f,.o, $(subst .f90,.o, $(MODULES)))))

# Reset suffixes that we understand
.SUFFIXES:
.SUFFIXES: .f90 .f .mod .o .F .F90

# Default Rules, the module rule should be executed first in most instances,
# this way the .mod file ends up always in the correct spot
%.mod : %.f90 ; $(FC) -c $< $(MODULE_FLAG)$(@D) $(ALL_INCLUDE) $(ALL_FFLAGS) -o $*.o
%.mod : %.f   ; $(FC) -c $< $(MODULE_FLAG)$(@D) $(ALL_INCLUDE) $(ALL_FFLAGS) -o $*.o
%.o : %.f90 ;   $(FC) -c $< 					$(ALL_INCLUDE) $(ALL_FFLAGS) -o $@
%.o : %.f ;     $(FC) -c $< 					$(ALL_INCLUDE) $(ALL_FFLAGS) -o $@

# debug:
#     @echo $(OBJECTS)
#     @echo $(MODULE_FILES)
#     @echo $(MODULE_PATHS)
#     @echo $(MODULE_OBJECTS)

$(EXECUTABLE): $(MODULE_FILES) $(MODULE_OBJECTS) $(OBJECTS) $(MAKEFILE_LIST) ;
	$(LINK) $(MODULE_OBJECTS) $(OBJECTS) $(ALL_INCLUDE) $(ALL_FFLAGS) -o $(EXECUTABLE)
	
clean:
	-rm -f $(EXECUTABLE) $(MODULE_FILES) $(MODULE_OBJECTS) $(OBJECTS)
	
debug:
	$(MAKE) new FFLAGS="$(DEBUG_FLAGS)"
	
opt:
	$(MAKE) new FFLAGS="$(OPT_FLAGS)"
    
new:
	$(MAKE) clean
	$(MAKE) $(EXECUTABLE)
