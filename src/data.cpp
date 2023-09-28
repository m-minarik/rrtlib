#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <iomanip>

#include "data.h"
#include "RAPID.H"
#include "types.h"

#include "matrix.h" // from libcip

bool off_to_rapid(const std::string &filename, RAPID_model *model, BoundingBox &bb, const Matrix *R_init, const Matrix *t_init)
{
	std::ifstream input(filename);
	if (!input.is_open())
	{
		std::cout << "Cannot open " << filename << "\n";
		return false;
	}

	double x, y, z;
	double x_, y_, z_;
	double v1, v2, v3;

	int vertex_count, face_count, edge_count;
	int face_id = 0;

	std::string line;
	int line_no = 1;

	// Begin building model
	model->BeginModel();
	struct Vertex
	{
		Vertex(double x, double y, double z)
		{
			coords[0] = x;
			coords[1] = y;
			coords[2] = z;
		}
		double coords[3];
	};

	std::vector<Vertex> vertices;


	while (std::getline(input, line))
	{
		std::stringstream stream(line);
		if (line_no == 1)
		{
			// First line - OFF
		}
		else if (line_no == 2)
		{
			// Second line - # of verts, faces, and edges
			stream >> vertex_count >> face_count >> edge_count;
			vertices.reserve(vertex_count);
		}
		else if (line_no - 2 <= vertex_count)
		{
			// Vertices - x y z
			if (R_init == nullptr)
			{
				stream >> x >> y >> z;
			}
			else
			{
				stream >> x_ >> y_ >> z_;
				x = (*R_init).val[0][0] * x_ + (*R_init).val[0][1] * y_ + (*R_init).val[0][2] * z_ + (*t_init).val[0][0];
				y = (*R_init).val[1][0] * x_ + (*R_init).val[1][1] * y_ + (*R_init).val[1][2] * z_ + (*t_init).val[0][1];
				z = (*R_init).val[2][0] * x_ + (*R_init).val[2][1] * y_ + (*R_init).val[2][2] * z_ + (*t_init).val[0][2];
			}

			vertices.emplace_back(x, y, z);

			if (x < bb.xmin)
				bb.xmin = x;
			if (x > bb.xmax)
				bb.xmax = x;
			if (y < bb.ymin)
				bb.ymin = y;
			if (y > bb.ymax)
				bb.ymax = y;
			if (z < bb.zmin)
				bb.zmin = z;
			if (z > bb.zmax)
				bb.zmax = z;
		}
		else
		{
			// Faces - 3 v1_id v2_id v3_id
			stream >> x >> v1 >> v2 >> v3;
			model->AddTri(vertices.at(int(v1)).coords,
						  vertices.at(int(v2)).coords,
						  vertices.at(int(v3)).coords,
						  face_id++);
		}
		++line_no;
	}

	// Finish building model
	model->EndModel();
	input.close();

	return true;
}
