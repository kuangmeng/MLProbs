#include <limits>
#include <cstdint>

#undef max
#undef min

// Generic sparse entry
template <class column_type, class value_type>
class SparseEntry
{
public:
	int getColumn() const	{ return first; }
	void setColumn(int c)	{ first = c; }

	float getValue() const	{ return second; }
	void setValue(float v)	{ second = v; }

protected:
	column_type first;
	value_type second;
};

// Fixed point specialisation
template <>
class SparseEntry<uint16_t, uint16_t>
{
public:
	int getColumn() const	{ return first; }
	void setColumn(int c)	{ first = c; }

	float getValue() const	{ return (float)second / std::numeric_limits<uint16_t>::max(); }
	void setValue(float v)	{ second = (uint16_t)(v * std::numeric_limits<uint16_t>::max()); }

protected:
	uint16_t first;
	uint16_t second;
};
