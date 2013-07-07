#ifndef __MANYCLAW_OCL_MANAGER_H
#define __MANYCLAW_OCL_MANAGER_H

#include <CL/cl.hpp>

class CLManager
{
public:
  CLManager(int argc, char **argv);

  void display();

  void addSourceFile(std::string filename);

  void addSource(std::string filename);

  void addArg(cl_mem_flags mem_flags, size_t size, void* host_ptr);

  void execute();

private:
  cl::Platform platforms;
  cl::Device devices;
  cl::Context context;
  cl::CommandQueue queue;
  cl::Program program;

  std::vector<cl::Buffer> args;
};

#endif  // __MANYCLAW_OCL_MANAGER_H
