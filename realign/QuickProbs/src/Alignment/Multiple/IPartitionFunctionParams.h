#pragma once

namespace quickprobs
{

class IPartitionFunctionParams
{
public:
	virtual ::size_t sizeInBytes() const = 0;
	virtual void* data() const = 0;
};

};