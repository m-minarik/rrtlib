#include <string>
#include <fstream>

#include "matrix.h"
#include "icpPointToPlane.h"

size_t load_object(const std::string &filename, double **object);

double getInitialTransformationICP(const std::string &model_path, const std::string &template_path, const std::string &correspondences_path,
                                   Matrix &R, Matrix &t)
{
	// Returns matrix R and vector t that transform the template onto the model

	  // allocate model and template memory
	double *M = nullptr;
	double *T = nullptr;

	std::cout << "Loading model from " << model_path << std::endl;
	size_t M_num = load_object(model_path, &M);
	std::cout << "Loading template from " << template_path << std::endl;
	size_t T_num = load_object(template_path, &T);
	std::cout << "Done\n";
	
	IcpPointToPlane icp(M, M_num, 3, "");
	double d;

	  // Load correspondences
	std::vector<std::pair<int, int> > correspondences;

	std::ifstream input(correspondences_path);
	if (!input.is_open() )
	{
		std::cout << "Cannot open " << correspondences_path << "\n";
		return 0;
	}

	int model_vert, template_vert;

	std::string line;

	while (std::getline(input, line) ) {
		std::stringstream stream (line);
        // Because the identify-object executable generates a file with
        // correspondences in the form
        // 'object' 'guiding_object'
        // which is 'template' 'model' in our case,
        // we need to swap the indices to be in the 'model' 'template' order
        // which is expected
		stream >> template_vert >> model_vert;
		correspondences.push_back( std::make_pair(model_vert, template_vert) );
	}

	  // Find transformation minimizing mutual distance between corresponding
	  // vertices
	d = icp.getInitialGuess(T, T_num, correspondences, R, t);

	  // run point-to-plane ICP (-1 = no outlier threshold)
	d = icp.fit(T, T_num, R, t, -1);

	  // free memory
	free(M);
	free(T);

	return d;
}

size_t load_object(const std::string &filename, double **object)
{
	std::ifstream input(filename);
	if (!input.is_open() )
	{
		std::cout << "Cannot open " << filename << "\n";
		return false;
	}

	double x, y, z;
	int vertex_count, face_count, edge_count;

	std::string line;
	int line_no = 1;
	size_t vertex_id = 0;

	while (std::getline(input, line) ) {
		std::stringstream stream (line);
		if (line_no == 1) {
			// First line - OFF
		} else if (line_no == 2) {
			  // Second line - # of verts, faces, and edges
			stream >> vertex_count >> face_count >> edge_count;
			*object = (double *) calloc( 3 * vertex_count, sizeof(double) );
		} else if (vertex_id < vertex_count) {
			  // Vertices - x y z
			stream >> x >> y >> z;
			(*object)[3 * vertex_id + 0] = x;
			(*object)[3 * vertex_id + 1] = y;
			(*object)[3 * vertex_id + 2] = z;
			++vertex_id;
		} else {
			// Faces - 3 v1_id v2_id v3_id
		}
		++line_no;
	}

	  // Finish building model
	input.close();

	return vertex_id;
}
