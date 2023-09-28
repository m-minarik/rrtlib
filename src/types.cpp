#include "types.h"
#include <string>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <cmath>
#include <iomanip>
#include <fstream>
using namespace std;

ostream &operator<<(ostream &os, const Point3 &r)
{
	os << '(' << r.x << ',' << r.y << ',' << r.z << ')';
	return os;
}

ostream &operator<<(ostream &os, const Point3_3 &r)
{
	os << '(' << r.x << ',' << r.y << ',' << r.z << ")\n(" << r.alpha << ',' << r.beta << ',' << r.gamma << ')';
	return os;
}

double get_random(const double from, const double to)
{
	return (double)(1.0 * rand() / (1.0 * RAND_MAX)) * (to - from) + from;
}

std::string align_to(std::string s, const unsigned int a, char c, bool ljust)
{
	size_t len = s.length();
	if (len > a)
		return s;

	if (ljust)
	{
		for (unsigned int j = 0; j < a - len; j++)
		{
			s = s + c;
		}
	}
	else
	{
		for (unsigned int j = 0; j < a - len; j++)
		{
			s = c + s;
		}
	}
	return s;
}

std::string align_to(double d, int a, bool ljust)
{
	std::ostringstream out;
	out.precision(3);
	out << std::fixed << d;
	return align_to(out.str(), a, ' ', ljust);
}
