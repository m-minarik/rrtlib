#include <iostream>
#include <cmath>

#include "types.h"
#include "object.h"
#include "data.h"

#include "matrix.h" // from libcip

Object::Object(const std::string &_filename, const Matrix *R_init, const Matrix *t_init) : filename{_filename}, bb{}
{
#ifdef DEBUG
	std::cout << "\n--- Loading object " << this->filename << " ---\n";
#endif

	this->model = new RAPID_model;

	bool successful = off_to_rapid(this->filename, this->model, this->bb, R_init, t_init);
	if (!successful)
	{
		std::cerr << "Failed loading RAPID model\n";
		exit(EXIT_FAILURE);
	}

#ifdef DEBUG
	std::cout << "Sucessfuly loaded RAPID model\n";
	std::cout << "Bounding box:\n";
	std::cout << "xmin: " << bb.xmin << "\t"
			  << "xmax: " << bb.xmax << "\n"
			  << "ymin: " << bb.ymin << "\t"
			  << "ymax: " << bb.ymax << "\n"
			  << "zmin: " << bb.zmin << "\t"
			  << "zmax: " << bb.zmax << "\n";
#endif

	// Initialize transformation matricess
	// omit this-> to save space
	R[0][0] = 1.0;
	R[0][1] = 0.0;
	R[0][2] = 0.0;
	R[1][0] = 0.0;
	R[1][1] = 1.0;
	R[1][2] = 0.0;
	R[2][0] = 0.0;
	R[2][1] = 0.0;
	R[2][2] = 1.0;

	T[0] = 0.0;
	T[1] = 0.0;
	T[2] = 0.0;
}

Object::~Object()
{
	delete this->model;
}

std::ostream &operator<<(std::ostream &os, const Object &obj)
{
	// Output current transformation matrix
	os << "| " << obj.R[0][0] << " " << obj.R[0][1] << " " << obj.R[0][2] << " | " << obj.T[0] << " |\n"
	   << "| " << obj.R[1][0] << " " << obj.R[1][1] << " " << obj.R[1][2] << " | " << obj.T[1] << " |\n"
	   << "| " << obj.R[2][0] << " " << obj.R[2][1] << " " << obj.R[2][2] << " | " << obj.T[2] << " |\n";
	return os;
}

void Object::transform(const Point3_3 &p)
{
	// Initialize transformation matrices
	// omit this-> to save space

	// http://lavalle.pl/planning/node104.html
	R[0][0] = cos(p.alpha) * cos(p.beta);
	R[0][1] = cos(p.alpha) * sin(p.beta) * sin(p.gamma) - sin(p.alpha) * cos(p.gamma);
	R[0][2] = cos(p.alpha) * sin(p.beta) * cos(p.gamma) + sin(p.alpha) * sin(p.gamma);

	R[1][0] = sin(p.alpha) * cos(p.beta);
	R[1][1] = sin(p.alpha) * sin(p.beta) * sin(p.gamma) + cos(p.alpha) * cos(p.gamma);
	R[1][2] = sin(p.alpha) * sin(p.beta) * cos(p.gamma) - cos(p.alpha) * sin(p.gamma);

	R[2][0] = -sin(p.beta);
	R[2][1] = cos(p.beta) * sin(p.gamma);
	R[2][2] = cos(p.beta) * cos(p.gamma);

	T[0] = p.x;
	T[1] = p.y;
	T[2] = p.z;
}
