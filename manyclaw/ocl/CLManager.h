#ifndef __MANYCLAW_OCL_MANAGER_H
#define __MANYCLAW_OCL_MANAGER_H

#ifdef HAVE_OCL

#define __CL_ENABLE_EXCEPTIONS

#include <CL/cl.hpp>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <utility>

class CLManager
{
public:

  CLManager();

  CLManager(int argc, char **argv);

  void display();

  void addSourceFile(std::string filename, std::string kernel_name);

  void addSource(std::string source, std::string kernel_name);

  void addArg(cl_mem_flags mem_flags, size_t size, void* host_ptr);

  void execute(int num_threads);

  void copyBackArg(int arg_number);

private:
  cl::Platform platform;
  std::vector<cl::Device> devices;
  cl::Context context;
  cl::CommandQueue queue;
  cl::Program program;
  cl::Kernel kernel;

  std::vector<cl::Buffer> arg_bufs;
  std::vector<size_t> arg_sizes;
  std::vector<void*> arg_ptrs;
};

#endif // HAVE_OCL

#endif  // __MANYCLAW_OCL_MANAGER_H
