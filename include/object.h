#ifndef __OBJECT__
#define __OBJECT__

#include <iostream>
#include <string>

#include "types.h"
#include "RAPID.H"

#include "matrix.h" // from libcip

struct Object
{
	Object(const std::string &_filename, const Matrix *R_init = nullptr, const Matrix *t_init = nullptr);
	~Object();

	void transform(const Point3_3 &p);

	friend std::ostream &operator<<(std::ostream &os, const Object &r);

	const std::string &filename;

	BoundingBox bb;

	double R[3][3];
	double T[3];
	RAPID_model *model;
};

#endif // __OBJECT__
