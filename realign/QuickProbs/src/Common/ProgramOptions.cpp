#include "ProgramOptions.h"

#include <list>
#include <algorithm>
#include <iomanip>

using namespace std;

bool ProgramOptions::parse(int argc, char *argv[])
{
	// copy parameters to temporary collection
	std::list<std::string> args(argc - 1);
	std::transform(argv + 1, argv + argc, args.begin(), [](char * a)->std::string { return std::string(a); });

	// parse normal options
	for (auto it = args.begin(); it != args.end();) {
		std::string param = *it;

		if (param[0] != '-') {
			++it;
			continue;
		}

		while (param[0] == '-') {
			param = param.substr(1);
		}

		if (options.find(param) != options.end()) {
			it = args.erase(it); // remove parsed option from collection
			auto option = options[param];

			// check if option is a switch
			auto switchOption = std::dynamic_pointer_cast<Switch>(option);
			if (switchOption) {
				switchOption->set(true);

			}
			else {
				const std::string& value = *it;
				if (option->parse(value)) {
					it = args.erase(it); // remove parsed value from collection
				}
				else {
					// unable to parse

				}
			}
		}
	}

	// everything left are positional options
	int posId = 0;

	for (auto it = args.begin(); it != args.end() && posId != positionalOptions.size(); ++posId, ++it) {
		std::string param = *it;
		positionalOptions[posId]->parse(param);
	}

	if (posId < positionalOptions.size()) {
		return false;
	}

	return true;
}


std::string ProgramOptions::toString(bool showDeveloperOptions) const
{
	std::ostringstream out;

	// print command line
	out << "Options:" << std::endl;
	for (const auto& o : options) {
		if (!o.second->getIsDeveloper() && shortNames.find(o.first) == shortNames.end()) {
			out << "\t" << std::setw(maxNameLength) <<  std::left << o.second->getName() << "\t" << o.second->getDescription() << std::endl;
		}
	}

	if (showDeveloperOptions) {
		out << "Developer options:" << std::endl;
		for (const auto& o : options) {
			if (o.second->getIsDeveloper() && shortNames.find(o.first) == shortNames.end()) {
				out << "\t" << std::setw(maxNameLength) << std::left << o.second->getName() << "\t" << o.second->getDescription() << std::endl;
			}
		}
	}

	return out.str();
}