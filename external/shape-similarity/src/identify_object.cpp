#include "GA.h"
#include "input_parser.h"

#include <string>
#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;

std::string LIBRARY = "";
std::string OBJECT = "";
std::string OUTPUT = "";
std::string TEMP = "temp";

bool parse_args(int argc, char **argv);

double evaluate_similarity(const std::string &a_file, int a_id, const std::string &b_file, int b_id);
void compute_initial_robust_map(const std::string &a_file, int a_id, const std::string &b_file, int b_id);
double compute_final_map(const std::string &a_file, int a_id, const std::string &b_file, int b_id);

int main(int argc, char **argv)
{
	// Parse command line arguments
	if (!parse_args(argc, argv))
	{
		return EXIT_FAILURE;
	}

	// Clear the temp folder
	system(("rm -r " + TEMP + "/* ; mkdir -p " + TEMP + "/robustMaps").c_str());

	double min_similarity = 0;
	std::string min_label = "";
	int min_index = 0;

	int index = 1;

	// Scans library for available objects, comparing the object similarity
	for (const auto &entry : fs::directory_iterator(LIBRARY))
	{
		// For each object directory, evaluate similarity
		std::string label = entry.path().filename();
		std::string library_object = entry.path().string() + "/object.off";

		double similarity = evaluate_similarity(OBJECT, 0, library_object, index);
		std::cout << label << ":\t" << similarity << "\n";

		if (similarity < min_similarity || index == 1)
		{
			min_similarity = similarity;
			min_label = label;
			min_index = index;
		}
		++index;
	}

	// Copy the result to OUTPUT and save other info
	// tail -n +2 "temp/0-1 ga+as.dat" "output/correspondences" - copy the correspondences without the first line, which contains metrics
	system(("tail -n +2 \"" + TEMP + "/0-" + std::to_string(min_index) + " ga+as.dat\" > \"" + OUTPUT + "/correspondences.txt\"").c_str());
	// ln -s "library/map/object/object.off" "output/guiding_object.off" - copy the guiding object symlink
	system(("ln -s \"" + LIBRARY + "/" + min_label + "/object.off\" \"" + OUTPUT + "/guiding_object.off\"").c_str());
	// ln -s "library/map/object" "output/guiding_paths"
	system(("ln -s \"" + LIBRARY + "/" + min_label + "\" \"" + OUTPUT + "/guiding_paths\"").c_str());

	std::ofstream output(OUTPUT + "/results.txt", std::ofstream::out);
	output << min_label << "\n"
		   << "\"" << LIBRARY + "/" + min_label << "\"\n"
		   << min_similarity << "\n";
	output.close();

	// Clear the TEMP folder
	system(("rm -r " + TEMP + "/*").c_str());
}

double evaluate_similarity(const std::string &a_file, int a_id, const std::string &b_file, int b_id)
{
	compute_initial_robust_map(a_file, a_id, b_file, b_id);
	return compute_final_map(a_file, a_id, b_file, b_id);
}

void compute_initial_robust_map(const std::string &a_file, int a_id, const std::string &b_file, int b_id)
{
	Mesh mesh_a{false, a_id, TEMP};
	Mesh mesh_b{true, b_id, TEMP};

	mesh_a.loadOff((char *)a_file.c_str());
	mesh_b.loadOff((char *)b_file.c_str());

	GA ga_initial(&mesh_a, &mesh_b, COMPUTE_ROBUST_LIST, 10, TEMP);
	ga_initial.go(false, true);
}

double compute_final_map(const std::string &a_file, int a_id, const std::string &b_file, int b_id)
{
	Mesh mesh_a{false, a_id, TEMP};
	Mesh mesh_b{true, b_id, TEMP};

	mesh_a.loadOff((char *)a_file.c_str());
	mesh_b.loadOff((char *)b_file.c_str());

	GA ga_final(&mesh_a, &mesh_b, GA_WITH_ROBUST_LIST, 10, TEMP);
	ga_final.go(false, true);

	return mesh_a.adaptiveSampling3(&mesh_b, true, TEMP);
}

bool parse_args(int argc, char **argv)
{
	InputParser input(argc, argv);

	if (input.cmdOptionExists("-h"))
	{
		std::cout << "--- Genetic Algorithms similarity evaluation ---\n\n"
				  << "Usage: \n\n"
				  << "-o"
				  << "\t"
				  << "OBJECT"
				  << "\t\t"
				  << "Object file (.off) location\n"
				  << "-l"
				  << "\t"
				  << "LIBRARY"
				  << "\t\t"
				  << "Library directory location\n"
				  << "-t"
				  << "\t"
				  << "TEMP"
				  << "\t\t"
				  << "Directory to which temporary files will be saved\n"
				  << "--out"
				  << "\t"
				  << "OUTPUT"
				  << "\t"
				  << "Directory to save the output to\n\n"
				  << "-h"
				  << "\t\t\t"
				  << "Show this help and exit"
				  << "\n";
		return false;
	}

	OBJECT = input.getCmdOption("-o");
	if (OBJECT.empty())
	{
		std::cout << "Path to an object file needs to be specified by the -o parameter\n";
		return false;
	}

	LIBRARY = input.getCmdOption("-l");
	if (LIBRARY.empty())
	{
		std::cout << "Library directory needs to be specified by the -l parameter\n";
		return false;
	}

	if (input.cmdOptionExists("-t"))
	{
		TEMP = input.getCmdOption("-t") + "/ga";
	}

	OUTPUT = input.getCmdOption("--out");
	if (OUTPUT.empty())
	{
		std::cout << "Output file needs to be specified by the --out parameter\n";
		return false;
	}

	return true;
}
