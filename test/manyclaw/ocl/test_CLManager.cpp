#include "gtest/gtest.h"

#include <manyclaw/manyclaw.h>

#include <limits.h>

const std::string vecAddSource =
    "__kernel void"
    "vadd(__global int * a, __global int * b, __global int *c)"
    "{"
    "    size_t i = get_global_id(0);"
    "    c[i] = a[i] + b[i];"
    "}";

TEST(CLManager, vecAdd) {
  CLManager manager(0, NULL);
  manager.addSource(vecAddSource);
  size_t buf_size = 3;
  int A[] = {0, 1, 2};
  int B[] = {0, 2, 3};
  int C[] = {0, 0, 0};
  manager.addArg(CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                 buf_size*sizeof(int), (void *) A);
  manager.addArg(CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                 buf_size*sizeof(int), (void *) B);
  manager.addArg(CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR,
                 buf_size*sizeof(int), (void *) C);
  manager.execute();
  manager.copyBackArg(3);
  EXPECT_EQ(0, C[0]);
  EXPECT_EQ(3, C[1]);
  EXPECT_EQ(5, C[2]);
}

