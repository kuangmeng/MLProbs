#pragma once
#include <map>
#include <queue>

template <class T, class U>
class FirstComparer
{
public:
	bool operator()(const std::pair<T, U*>& a, const std::pair<T, U*>& b) {
		return a.first < b.first;
	}
};

template <class T>
class MemoryPool
{
public:
	MemoryPool(int count, size_t size) {
		for (int i = 0; i < count; ++i) {
			T* a = new T[size];
			areas.insert(std::pair<int,T*>(std::numeric_limits<int>::max(), a));
		}
	}

	~MemoryPool() {
		for (auto& a : areas) {
			delete [] a.second;
		}
	}
	
	bool isAssigned(void* owner) {
		return owners2Areas.find(owner) != owners2Areas.end();
	}

	T* assign(void *owner, int priority) {
		for (auto& a : areas) {
			if (priority <= a.first) {
				T* area = a.second; 
				owners2Areas[owner] = area;
				areas2Owners[area] = owner;
				return area;
			}
		}
	} 
	
	void* get(void* owner) {
		return owners2Areas.at(owner);
	}

protected:
	std::multimap<int,T*> areas;
	std::map<void*, T*> owners2Areas;
	std::map<T*, void*> areas2Owners; 
};