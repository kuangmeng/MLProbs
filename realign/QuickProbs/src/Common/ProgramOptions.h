#pragma once
#undef max

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <memory>
#include <map>
#include <algorithm>
#include <set>

class AbstractOption {
public:
	std::string getName() const { return name; }
	std::string getDescription() const { return desc; }
	bool getIsDeveloper() const { return isDeveloper; }

	AbstractOption(const std::string& name, const std::string& desc, bool isDeveloper) : name(name), desc(desc), isDeveloper(isDeveloper) {}

	virtual bool parse(const std::string& valstr) = 0;

protected:
	const std::string name;
	const std::string desc;
	const bool isDeveloper;
};

template <class T>
class Option : public AbstractOption {
public:
	bool getIsSet() const { return isSet; }
	
	Option(const std::string& name, const std::string& desc, bool isDeveloper)
		: AbstractOption(name, desc, isDeveloper), isSet(false) {}

	Option(const std::string& name, const std::string& desc, bool isDeveloper, const T def)
		: AbstractOption(name, desc, isDeveloper), value(def), isSet(true) {}

	T get() const { return value; }
	void set(T v) { this->value = v; isSet = true;  }

	virtual bool parse(const std::string& valstr) {
		std::istringstream iss(valstr);
		bool ok = static_cast<bool>(iss >> value);
		isSet |= ok;
		return ok;
	}

protected:
	T value;
	bool isSet ;
};

class Switch : public Option<bool>
{
public:
	Switch(const std::string& name, const std::string& desc, bool isDeveloper) : Option(name, desc, isDeveloper, false) {}
};


class ProgramOptions
{
public:
	ProgramOptions(const std::string & executable) : executable(executable), maxNameLength(0)
	{

	}

	bool exists(const std::string & name) const
	{
		return options.find(name) != options.end();
	}

	/// <summary>
	/// Adds option with given name and description.
	/// </summary>
	/// <param name ="name">Option name.</parameter>
	/// <param name ="desc">Option description.</parameter>
	template <typename Type>
	void add(const std::string & name, const std::string & desc, Type def, bool developer)
	{  
		size_t sepPos = name.find(',');
		std::string mainName = (sepPos == std::string::npos) ? name : name.substr(0,sepPos);
		std::string shortName = (sepPos == std::string::npos) ? "" : name.substr(sepPos + 1);

		auto option = std::shared_ptr<AbstractOption>(new Option<Type>(name, desc, developer, def));
		options[mainName] = option;
		if (shortName.length() > 0) {
			options[shortName] = option;
			shortNames.insert(shortName);
		}
		maxNameLength = std::max(maxNameLength, name.length());
	}

	template <typename Type>
	void add(const std::string & name, const std::string & desc, bool developer)
	{
		size_t sepPos = name.find(',');
		std::string mainName = (sepPos == std::string::npos) ? name : name.substr(0, sepPos);
		std::string shortName = (sepPos == std::string::npos) ? "" : name.substr(sepPos + 1);

		auto option = std::shared_ptr<AbstractOption>(new Option<Type>(name, desc, developer));
		options[mainName] = option;
		if (shortName.length() > 0) {
			options[shortName] = option;
			shortNames.insert(shortName);
		}

		maxNameLength = std::max(maxNameLength, name.length());
	}

	
	template <typename Type>
	void addPositional(const std::string & name, const std::string & desc)
	{
		auto option = std::shared_ptr<AbstractOption>(new Option<Type>(name, desc, false));
		positionalOptions.push_back(option);

		maxNameLength = std::max(maxNameLength, name.length());
	}

	void addSwitch(const std::string & name, const std::string & desc, bool developer)
	{
		size_t sepPos = name.find(',');
		std::string mainName = (sepPos == std::string::npos) ? name : name.substr(0, sepPos);
		std::string shortName = (sepPos == std::string::npos) ? "" : name.substr(sepPos + 1);

		auto option = std::shared_ptr<AbstractOption>(new Switch(name, desc, developer));
		options[mainName] = option;
		if (shortName.length() > 0) {
			options[shortName] = option;
			shortNames.insert(shortName);
		}

		maxNameLength = std::max(maxNameLength, name.length());
	}


	bool parse(int argc, char *argv[]);

	
	template <typename T>
	bool get(const std::string & name, T &result)
	{
		auto it = options.find(name);
		if (it != options.end()) {
			auto option = std::dynamic_pointer_cast<Option<T>>(it->second);
			if (option->getIsSet()) {
				result = option->get();
				return true;
			}
		}

		auto it2 = std::find_if(positionalOptions.begin(), positionalOptions.end(), [&name](std::shared_ptr<AbstractOption>& o) { return o->getName() == name; });
		if (it2 != positionalOptions.end()) {
			auto option = std::dynamic_pointer_cast<Option<T>>(*it2);
			if (option->getIsSet()) {
				result = option->get();
				return true;
			}
		}
			
		return false;
	}

	template <typename T>
	T get(const std::string & name)
	{
		T var;
		if (!get(name, var)) {
			throw std::runtime_error("Unable to read parameter: " + name);
		}

		return var;
	}

	std::string toString(bool showDeveloperOptions) const;

private:
	std::map<std::string, std::shared_ptr<AbstractOption>> options;
	std::set<std::string> shortNames;
	std::vector<std::shared_ptr<AbstractOption>> positionalOptions;

	size_t maxNameLength;

	std::string executable;
};