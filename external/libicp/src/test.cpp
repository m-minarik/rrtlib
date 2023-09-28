#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <utility>
#include "icpPointToPlane.h"

const std::string exec(const char *cmd);
size_t load_object(const std::string &filename, double **object, const Matrix &R, const Matrix &t);

// main `model_path` `template_path` `angle` `t` `correspondences`
int main (int argc, char **argv)
{
	  // Generate random rotation and translation
	Matrix R = Matrix::rotFromAngle( std::atof(argv[3]) );
	Matrix t(3, 1);

	t.val[0][0] = 3*std::atof(argv[4]);
	t.val[0][1] = std::atof(argv[4]);
	t.val[0][2] = -std::atof(argv[4]);

	  // allocate model and template memory
	double *M = nullptr;
	double *T = nullptr;

	std::cout << R << "\n";
	std::cout << t << "\n";

	size_t M_num = load_object(argv[1], &M, R, t);
	std::cout << "Loaded model " << argv[1] << " with " << M_num << " points ..." << "\n";

	  // start with identity as initial transformation
	  // in practice you might want to use some kind of prediction here
	R.eye();
	t.zero();

	size_t T_num = load_object(argv[2], &T, R, t);
	std::cout << "Loaded template " << argv[2] << " with " << T_num << " points ..." << "\n";


	  // Determine save folder name and create it
	  // int run_id = std::stoi( exec("ls ../tests/ | wc -l") );
	  // std::string save_dir{"../tests/" + std::to_string(run_id + 1) + "/"};
	std::string save_dir{""};
	//
	// int res = std::system( ("mkdir -p " + save_dir).c_str() );
	// if (res != 0) {
	// 	std::cerr << "Couldn't create " << save_dir << "\n";
	// 	exit(EXIT_FAILURE);
	// }

	IcpPointToPlane icp(M, M_num, 3, save_dir);
	double d;
	  // Load correspondences
	if (argc == 6) {
		std::vector<std::pair<int, int> > correspondences;

		std::ifstream input(argv[5]);
		if (!input.is_open() )
		{
			std::cout << "Cannot open " << argv[5] << "\n";
			return 0;
		}

		int model_vert, template_vert;

		std::string line;

		while (std::getline(input, line) ) {
			std::stringstream stream (line);
			stream >> model_vert >> template_vert;
			correspondences.push_back( std::make_pair(model_vert, template_vert) );
		}

		d = icp.getInitialGuess(T, T_num, correspondences, R, t);

		std::cout << "Initial guess:" << "\n";
		std::cout << "R:" << "\n" << R << "\n" << "\n";
		std::cout << "t:" << "\n" << t << "\n" << "\n";
	}

	  // run point-to-plane ICP (-1 = no outlier threshold)
	// std::cout << "Running ICP (point-to-plane, no outliers)" << "\n";
	d = icp.fit(T, T_num, R, t, -1);

	// std::ofstream cmd(save_dir + "cmd.sh");
	// std::ofstream results(save_dir + "results.txt");
	//
	// for (int i = 0; i < argc; ++i){
	// 	cmd << argv[i] << " ";
	// }
	//
	// results << R << "\n\n" << t;

	  // results
	std::cout << "Transformation results:" << "\n";
	std::cout << "R:" << "\n" << R << "\n" << "\n";
	std::cout << "t:" << "\n" << t << "\n" << "\n";

	std::cout << "Final difference=(" << d << ")\n";
	// cmd.close();
	// results.close();

	  // free memory
	free(M);
	free(T);

	  // success
	return 0;
}

size_t load_object(const std::string &filename, double **object, const Matrix &R, const Matrix &t)
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
			(*object)[3 * vertex_id + 0] = R.val[0][0] * x + R.val[0][1] * y + R.val[0][2] * z + t.val[0][0];
			(*object)[3 * vertex_id + 1] = R.val[1][0] * x + R.val[1][1] * y + R.val[1][2] * z + t.val[0][1];
			(*object)[3 * vertex_id + 2] = R.val[2][0] * x + R.val[2][1] * y + R.val[2][2] * z + t.val[0][2];
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

// https://stackoverflow.com/questions/478898/how-do-i-execute-a-command-and-get-the-output-of-the-command-within-c-using-po
const std::string exec(const char *cmd)
{
	std::array<char, 128> buffer;
	std::string result;
	std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
	if (!pipe) {
		std::cerr << "popen() failed!";
		exit(EXIT_FAILURE);
	}
	while (fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr) {
		result += buffer.data();
	}
	return result;
}
