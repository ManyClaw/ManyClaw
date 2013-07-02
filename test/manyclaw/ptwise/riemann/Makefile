# For simplicity and to avoid depending on Google Test's
all : execute_test

MANYCLAW_DIR=../../../..
include $(MANYCLAW_DIR)/variables.inc

TEST_SRC = $(wildcard *.cpp)
TEST_OBJ =  $(patsubst %.cpp,%.o,$(TEST_SRC))

include $(MANYCLAW_DIR)/test/test_rules.inc
