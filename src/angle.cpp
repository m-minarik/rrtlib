#include "angle.h"
#include <iostream>

Angle operator+(const Angle &lhs, const Angle &rhs)
{
	return Angle(lhs.value + rhs.value);
}

Angle operator+(const Angle &lhs, const double &rhs)
{
	return Angle(lhs.value + rhs);
}

double operator-(const Angle &lhs, const Angle &rhs)
{
	return wrap_to_pmPI(lhs.value - rhs.value);
}

Angle operator-(const Angle &lhs, const double &rhs)
{
	return Angle(lhs.value - rhs);
}
