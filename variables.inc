# Library settings
NAME ?= manyclaw
INSTALL_PATH ?= lib
MAJOR_VERSION ?= 0
MINOR_VERSION ?= 1.0


CXX = g++
LINK ?= $(CXX)
ifeq ($(CXX), g++)
	OMPCFLAGS = -fopenmp
else ifeq ($(CXX), icc)
	OMPCFLAGS = -openmp
endif


CXXFLAGS = -I./ -O0 -g -Wall $(OMPCFLAGS)
LDFLAGS = 


ifdef TBB_LIB
      CXXFLAGS += -I$(TBB_INC)
      TBB_LDFLAGS = -L$(TBB_LIB)  
endif
TBB_LDFLAGS += -ltbb
LDFLAGS += $(TBB_LD_FLAGS)


# Dynamic library flags (archetecture dependent), note the platforms here are
# found by sys.platform in ptyhon, check
UNAME = $(shell uname)
ifeq ($(UNAME), Linux)
    LIB_FULL_NAME = lib$(NAME).so.$(MAJOR_VERSION).$(MINOR_VERSION)
    LIB_MAJOR_NAME = lib$(NAME).so.$(MAJOR_VERSION)
    LIB_SHORT_NAME = lib$(NAME).so
    LIB_INSTALL_NAME = -Wl,-soname,$(LIB_MAJOR_NAME)
    VERSION_FLAGS =
    
    # Compiler specific flags
    ifeq ($(CXX), g++)
        CXXFLAGS += -fPIC
        LDFLAGS += -shared
    endif
else 
    ifeq ($(UNAME), Darwin)
        # Darwin, has a few more options and calls its libraries something
        # different but is taken care of here.  Note that we have to use the 
        # flag -flat_namespace in order for it to act like linux, see the ld 
        # man page for more information on the two-level namespace OS X uses
        LIB_FULL_NAME = lib$(NAME).$(MAJOR_VERSION).$(MINOR_VERSION).dylib
        LIB_MAJOR_NAME = lib$(NAME).$(MAJOR_VERSION).dylib
        LIB_SHORT_NAME = lib$(NAME).dylib
        LIB_INSTALL_NAME = -install_name $(LIB_MAJOR_NAME)
        VERSION_FLAGS = -compatibility_version $(MAJOR_VERSION) -current_version $(MAJOR_VERSION).$(MINOR_VERSION)
    
        # Compiler specific flags
        ifeq ($(CXX), g++)
            CXXFLAGS += -fPIC
            LDFLAGS += -dynamiclib -flat_namespace -undefined suppress
        else 
            ifeq ($(CXX), g++)
            	CXXFLAGS += -module $(MOD_PATH)
                LDFLAGS += -dynamiclib -flat_namespace -undefined suppress
            endif
        endif
    endif
endif