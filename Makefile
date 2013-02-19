include variables.inc

.PHONY.: all, build, lib, test, benchmark, clean, clobber

all: lib install test benchmark

build:
	$(MAKE) -C manyclaw/common
	$(MAKE) -C manyclaw/ptwise/steppers/advection
	$(MAKE) -C manyclaw/ptwise/steppers/acoustics_const
	$(MAKE) -C manyclaw/ptwise/steppers/acoustics_var
	$(MAKE) -C manyclaw/ptwise/steppers/euler

lib: $(LIB_FULL_NAME)

$(LIB_FULL_NAME): build
	$(LINK) $(LIB_INSTALL_NAME) $(VERSION_FLAGS) $(shell find manyclaw -name "*.o") -o $(LIB_FULL_NAME) $(INCLUDE) $(LDFLAGS)

install: $(INSTALL_PATH)/$(LIB_FULL_NAME)

$(INSTALL_PATH)/$(LIB_FULL_NAME): $(LIB_FULL_NAME)
	@if [ ! -d ${INSTALL_PATH} ]; then \
	  mkdir $(INSTALL_PATH); \
	fi
	-cp $(LIB_FULL_NAME) $(INSTALL_PATH)
	@if [ ! -e $(INSTALL_PATH)/$(LIB_MAJOR_NAME) ]; then \
	  ln -s $(LIB_FULL_NAME) $(INSTALL_PATH)/$(LIB_MAJOR_NAME); \
	fi
	@if [ ! -e $(INSTALL_PATH)/${LIB_SHORT_NAME} ]; then \
	  ln -s $(LIB_FULL_NAME) $(INSTALL_PATH)/$(LIB_SHORT_NAME); \
	fi

test: lib
	$(MAKE) -C test

benchmark: lib
	$(MAKE) -C benchmark

clean:
	$(MAKE) -C manyclaw/common clean
	$(MAKE) -C manyclaw/ptwise/steppers/advection clean
	$(MAKE) -C manyclaw/ptwise/steppers/acoustics_const clean
	$(MAKE) -C manyclaw/ptwise/steppers/acoustics_var clean
	$(MAKE) -C manyclaw/ptwise/steppers/euler clean
	$(MAKE) -C test clean
	$(MAKE) -C benchmark clean
	@echo $(LIB_FULL_NAME) $(INSTALL_PATH)
	-rm -rf $(INSTALL_PATH)/$(LIB_FULL_NAME)

clobber: clean
	-rm -rf $(INSTALL_PATH)/$(LIB_MAJOR_NAME)
	-rm -rf $(INSTALL_PATH)/$(LIB_SHORT_NAME)