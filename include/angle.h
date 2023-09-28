#ifndef __ANGLE__
#define __ANGLE__

#define _PI 3.141592653589793
#define _2PI 6.283185307179586
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

inline double wrap_to_PI(const double _value)
{
	double value = _value;

	while (value >= _PI)
	{
		value -= _PI;
	}

	while (value < 0)
	{
		value += _PI;
	}
	return value;
}

inline double wrap_to_2PI(const double _value)
{
	double value = _value;

	while (value >= _2PI)
	{
		value -= _2PI;
	}

	while (value < 0)
	{
		value += _2PI;
	}
	return value;
}

inline double wrap_to_pmPI(const double _value)
{
	double value = _value;

	while (value >= _PI)
	{
		value -= _2PI;
	}

	while (value < -_PI)
	{
		value += _2PI;
	}
	return value;
}

class Angle
{
public:
	Angle(const double &_value) : value{wrap_to_2PI(_value)} {};

	Angle &operator=(const Angle &other)
	{
		this->value = other.value;
		return *this;
	}

	Angle &operator=(const double &value)
	{
		this->value = wrap_to_2PI(value);
		return *this;
	}

	Angle &operator+=(const double &value)
	{
		this->value = wrap_to_2PI(this->value + value);
		return *this;
	}

	operator double() const
	{
		return this->value;
	}

	double value;
};

Angle operator+(const Angle &lhs, const Angle &rhs);
Angle operator+(const Angle &lhs, const double &rhs);
double operator-(const Angle &lhs, const Angle &rhs);
Angle operator-(const Angle &lhs, const double &rhs);

#endif
