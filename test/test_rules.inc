TEST_EXE = tester

TOP_DIR = $(MANYCLAW_DIR)/test
GTEST_DIR = $(MANYCLAW_DIR)/contrib

CPPFLAGS += -I$(GTEST_DIR) -I$(MANYCLAW_DIR)
LDFLAGS = -L$(MANYCLAW_DIR)/lib -lmanyclaw $(TBB_LDFLAGS)

# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/gtest/gtest.h 

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/gtest/gtest-all.cc $(GTEST_HEADERS)

# Builds gtest.a and gtest_main.a.o
# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
$(TOP_DIR)/gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/gtest/gtest-all.cc

$(TOP_DIR)/gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/gtest/gtest_main.cc

$(TOP_DIR)/gtest.a : $(TOP_DIR)/gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

$(TOP_DIR)/gtest_main.a : $(TOP_DIR)/gtest-all.o $(TOP_DIR)/gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

tester : $(TEST_OBJ) $(TOP_DIR)/gtest_main.a
	echo $(TEST_OBJ)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@ -lpthread $(LDFLAGS)

execute_test : tester
	./tester

clean :
	rm -f $(TEST_EXE) $(TEST_OBJ) $(TOP_DIR)/gtest.a $(TOP_DIR)/gtest_main.a $(TOP_DIR)/gtest_main.o *.o  $(call rwildcard,manyclaw,tester)

.PHONY: all execute_test clean clobber
