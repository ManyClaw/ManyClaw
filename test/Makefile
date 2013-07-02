MANYCLAW_DIR=..
include $(MANYCLAW_DIR)/variables.inc
GTEST_DIR = $(MANYCLAW_DIR)/contrib

CPPFLAGS += -I$(GTEST_DIR) -I$(MANYCLAW_DIR)
CXXFLAGS += -g -Wall -Wextra
LDFLAGS = -L$(MANYCLAW_DIR)/lib -lmanyclaw $(TBB_LDFLAGS)

# Make does not offer a recursive wildcard function, so here's one:
rwildcard=$(wildcard $1$2) $(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2))

TEST_SRC = $(call rwildcard,manyclaw,*.cpp)
TEST_OBJ =  $(patsubst %.cpp,%.o,$(TEST_SRC))
TEST_EXE = test_all

# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/gtest/gtest.h 

.PHONY = all execute_test clean clobber

all : execute_test

# Builds gtest.a and gtest_main.a.

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/gtest/gtest-all.cc $(GTEST_HEADERS)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/gtest/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/gtest/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

test_all : $(TEST_OBJ) gtest_main.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@ -lpthread $(LDFLAGS)

execute_test : test_all
	./test_all

clean :
	rm -f $(TEST_EXE) $(TEST_OBJ) gtest.a gtest_main.a *.o

clobber : clean
	rm -f $(TEST_EXE)
