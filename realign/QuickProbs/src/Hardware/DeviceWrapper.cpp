#include "DeviceWrapper.h"

#include "Common/Log.h"

clex::DeviceWrapper::DeviceWrapper(const cl::Device& ref, const cl::Context& context,  bool useKernelProfiling)
{
	LOG_DEBUG << "Creating device info...";
	this->device = ref;
	this->info = std::make_shared<clex::DeviceInfo>(ref);
	LOG_DEBUG << "done!" << std::endl;

	int props = 0;
	if (useKernelProfiling) {
		props |= CL_QUEUE_PROFILING_ENABLE;
	}

	// create five command queues
	LOG_DEBUG << "Creating queues...";
	int code;
	for (int i = 0; i < 10; ++i) {
		queues.push_back(std::make_shared<cl::CommandQueue>(context, device, props, &code));
		//clCall(code);
	}
	LOG_DEBUG << "done!" << std::endl;

	mainQueue = queues[queues.size() - 1];
}
