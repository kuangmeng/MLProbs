#include "StatisticsProvider.h"
#include <fstream>

std::ostream& operator<<(std::ostream& os, IStats& stats) {
	os << stats.toString();
	return os;
}

std::string StatisticsProvider::printStats()
{
	std::ostringstream oss;

	for (auto s = statistics.begin(); s != statistics.end(); s++) {
		oss << s->first << ": " << s->second->toString() << std::endl;
	}

	return oss.str();
}

void StatisticsProvider::saveStats(std::string file)
{
	std::ofstream out;
	out.open(file);
	
	out << "[stats]" << std::endl;

	for (auto s = statistics.begin(); s != statistics.end(); s++) {
		out << s->first << "=" << s->second->toString() << std::endl;
	}

	out.close();
}

void StatisticsProvider::joinStats(const StatisticsProvider& other)
{
	for (auto s = other.statistics.begin(); s != other.statistics.end(); ++s) {
		this->statistics.insert(*s);
	}
}
