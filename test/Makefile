all : execute_test

MANYCLAW_DIR=..
include $(MANYCLAW_DIR)/variables.inc

TEST_SRC = $(call rwildcard,manyclaw,*.cpp)
TEST_OBJ =  $(patsubst %.cpp,%.o,$(TEST_SRC))

include $(MANYCLAW_DIR)/test/test_rules.inc
