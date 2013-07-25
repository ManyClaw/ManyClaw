#include "CLManager.h"

#ifdef HAVE_OCL

CLManager::CLManager()
{
  try {
    std::vector<cl::Platform> platformList;
    cl::Platform::get(&platformList);
    platform = platformList[0];
    cl_context_properties cprops[] = {CL_CONTEXT_PLATFORM,
        (cl_context_properties) (platformList[0])(),
        0};
    context = cl::Context(CL_DEVICE_TYPE_CPU, cprops);
    devices = context.getInfo<CL_CONTEXT_DEVICES>();
    queue = cl::CommandQueue(context, devices[0], 0);
  } catch (cl::Error err) {
    std::cerr << "ERROR: " << err.what() << "(" << err.err() << ")" << std::endl;
    throw err;
  }
}

void CLManager::display()
{}

void CLManager::addSourceFile(std::string filename, std::string kernel_name)
{
  std::ifstream file(filename.c_str());
  std::string prog(std::istreambuf_iterator<char>(file),
                   (std::istreambuf_iterator<char>()));
  addSource(prog, kernel_name);
}

void CLManager::addSource(std::string program_src, std::string kernel_name)
{
  try {
    cl::Program::Sources sources(1, std::make_pair(program_src.c_str(), 0));
    program = cl::Program(context, sources);

    program.build(devices);
    kernel = cl::Kernel(program, kernel_name.c_str());
  } catch (cl::Error err) {
    std::cerr << "ERROR: " << err.what() << "(" << err.err() << ")" << std::endl;
    if (err.err() == CL_BUILD_PROGRAM_FAILURE) {
        std::string buildLog;
        program.getBuildInfo(devices[0], CL_PROGRAM_BUILD_LOG, &buildLog);
        std::cerr << "Error in kernel: " << std::endl
            << buildLog << std::endl;
    }
    throw err;
  }

}

void CLManager::addArg(cl_mem_flags mem_flags, size_t size, void* host_ptr)
{
  try {
      cl::Buffer arg_buf = cl::Buffer(context, mem_flags, size, host_ptr);
      arg_bufs.push_back(arg_buf);
      arg_sizes.push_back(size);
      arg_ptrs.push_back(host_ptr);
  } catch (cl::Error err) {
    std::cerr << "ERROR: " << err.what() << "(" << err.err() << ")" << std::endl;
    throw err;
  }
}

void CLManager::execute(int num_threads)
{
  try {
    cl::Event event;
    for (size_t i=0; i < arg_bufs.size(); ++i)
      kernel.setArg(i, arg_bufs[i]);
    queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(num_threads),
        cl::NullRange, NULL, &event);
    event.wait();
  } catch (cl::Error err) {
    std::cerr << "ERROR: " << err.what() << "(" << err.err() << ")" << std::endl;
    throw err;
  }
}

void CLManager::copyBackArg(int arg_number)
{
  try {
      queue.enqueueReadBuffer(arg_bufs[arg_number],
          CL_TRUE, //block
          0,
          arg_sizes[arg_number],
          arg_ptrs[arg_number]
      );
  } catch (cl::Error err) {
    std::cerr << "ERROR: " << err.what() << "(" << err.err() << ")" << std::endl;
    throw err;
  }
}

#endif /* HAVE_OCL */
