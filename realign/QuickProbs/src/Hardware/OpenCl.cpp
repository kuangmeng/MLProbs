#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "OpenCl.h"
#include "OpenCLError.h"
#include "PlatformInfo.h"
#include "DeviceInfo.h"

#include "../Common/dbgnew.h"
#include "../Common/Timer.h"
#include "../Common/MemoryTools.h"
#include "../Common/Log.h"


cl_int clCall(cl_int code)
{
	std::string msg = clex::OpenCLError::name(code);

	if(code != CL_SUCCESS) {
		std::cout << "Error: " << msg << std::endl;
	}

	return code;
}

double clex::OpenCL::profileTimeMsec(cl::Event & profilingEvent)
{
	cl_ulong start = profilingEvent.getProfilingInfo<CL_PROFILING_COMMAND_START>();
	cl_ulong end = profilingEvent.getProfilingInfo<CL_PROFILING_COMMAND_END>();
	double timeMsec = (double)(end - start) / 1000000.0;
	return timeMsec;
}

std::string clex::OpenCL::listDevices(int deviceType)
{
	std::ostringstream oss;
	std::vector<cl::Platform> platforms;
	cl::Platform::get(&platforms);

	for (int i = 0; i < platforms.size(); i++) {
		oss << "Platform " << i << ": " << PlatformInfo(platforms[i]).getName() << std::endl;
		std::vector<cl::Device> devices;
		platforms[i].getDevices(deviceType, &devices);

		for (int j = 0; j < devices.size(); j++) {
			oss << "\tDevice " << j << ": " << DeviceInfo(devices[j]).getName() << std::endl;
		}
	}

	return oss.str();
}

clex::OpenCL::OpenCL(int deviceType, int platformNum, int deviceNum, bool kernelProfiling) :
	deviceType(deviceType)
{
	LOG_DEBUG << "Getting platforms...";
	int code1, code2, code3;
	std::vector<cl::Platform> platforms;
	code1 = clCall(cl::Platform::get(&platforms));
	LOG_DEBUG << "done!" << std::endl;
	for (auto& p : platforms) {
		//LOG_DEBUG << "P = " << p() << std::endl;
		PlatformInfo pi(p);
		LOG_DEBUG << pi.getName() << std::endl;
	}

	LOG_DEBUG << "Getting devices...";
	std::vector<cl::Device> tempDevices;
	code2 = clCall(platforms[platformNum].getDevices(deviceType, &tempDevices));
	LOG_DEBUG << "done!" << std::endl;
	for (auto& d : tempDevices) {
		//LOG_DEBUG << "D = " << d() << std::endl;
		DeviceInfo di(d);
		LOG_DEBUG << di.getName() << std::endl;
	}
	
	// fixme:
	std::vector<cl::Device> finalDevices(1, tempDevices[deviceNum]);
	
	LOG_DEBUG << "Creating context...";
	context = std::shared_ptr<cl::Context>(new cl::Context(finalDevices, NULL, &OpenCL::errorCallback, NULL, &code3));
	clCall(code3);
	LOG_DEBUG << "done!" << std::endl;

	LOG_DEBUG << "Creating device wrapper...";
	this->devices.push_back(std::shared_ptr<DeviceWrapper>(new DeviceWrapper(finalDevices[0], *context, kernelProfiling)));
	this->mainDevice = devices[0];
	LOG_DEBUG << "done!" << std::endl;

	if (code1 != CL_SUCCESS || code2 != CL_SUCCESS || code3 != CL_SUCCESS) {
		throw std::runtime_error("OpenCL context initialisation error!");
	}
}


double clex::OpenCL::profileKernel(
		const cl::Kernel& kernel,
		const cl::NDRange& offset,
		const cl::NDRange& global,
		const cl::NDRange& local)
{
	cl::Event finishEvent;

	clCall(mainDevice->mainQueue->enqueueNDRangeKernel(kernel, offset, global, local, NULL, &finishEvent));

	cl_ulong timeStart, timeEnd;
	double timeMsec;
	finishEvent.wait();
	finishEvent.getProfilingInfo(CL_PROFILING_COMMAND_START, &timeStart);
	finishEvent.getProfilingInfo(CL_PROFILING_COMMAND_END, &timeEnd);
	timeMsec = (timeEnd - timeStart) / 1000000.0;

	return timeMsec;
}
