#pragma once
#include <string>
#include <map>
#include <sstream>
#include <cstdint>
#include <memory>

#define GET_STATS

#ifdef GET_STATS
	#define STATS_WRITE(key,value) writeStats(key, value)
	#define STATS_ADD(key,value) addStats(key, value)
#else
	#define STATS_WRITE(key,value) 
	#define STATS_ADD(key,value) 
#endif


class IStats {
public:
	virtual ~IStats() {}
	virtual std::string toString() = 0;
	virtual void add(const IStats& other) = 0;

	virtual std::shared_ptr<IStats> clone() = 0;
};

template <class T>
class Stats : public IStats {
public:
	typedef T value_type;
	

	Stats(T value) : value(value) {}
	
	std::shared_ptr<IStats> clone() {
		auto copy = std::make_shared<Stats<T>>(this->value);
		return copy;
	}
	
	virtual void add(T other) {
		value += other;
	}

	virtual void add(const IStats& other) {
		auto casted = dynamic_cast<const Stats<T>&>(other);
		value += casted.value;
	}

	virtual std::string toString() { return std::to_string(value); }
protected:
	T value;
}; 

std::ostream& operator<<(std::ostream& os, IStats& stats);

class StatisticsProvider
{
public:
	
	virtual ~StatisticsProvider() {}

	template<class T>
	void writeStats(std::string key, T value) {
		statistics[key] = std::make_shared<Stats<T>>(value);
	}

	template <class T>
	void addStats(std::string key, T other) {
		auto casted = std::dynamic_pointer_cast<Stats<T>>(statistics[key]);
		casted->add(other);
	}

	void addStats(const StatisticsProvider& other) {
		for (const auto &s : other.statistics) {
			if (this->statistics.find(s.first) == this->statistics.end()) {
				this->statistics[s.first] = s.second->clone();
			} else {
				this->statistics[s.first]->add(*s.second);
			}
		}
	}

	void joinStats(const StatisticsProvider& other);

	void clearStats() { statistics.clear(); }

	std::string printStats();

	void saveStats(std::string file);

protected:
	std::map<std::string, std::shared_ptr<IStats>> statistics;


};