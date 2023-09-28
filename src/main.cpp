#include <iostream>
#include <memory>
#include <cstdlib> // std::system
#include <filesystem>
#include <cmath>

#include <ompl/tools/benchmark/Benchmark.h>
#include <ompl/geometric/planners/rrt/RRTLIB.h>
#include <ompl/base/spaces/SE3StateSpace.h>

#include "validity_checker.h"
#include "icpWrapper.h"
#include "types.h"
#include "input_parser.h"
#include "ompl_utils.h"

namespace fs = std::filesystem;

// Default parameter values - can be overwritten by command line arguments
double PLANNER_RANGE = 0.1;
double GOAL_BIAS = 0.05;
double PATH_BIAS = 0.80;
double GUIDING_RADIUS = 0.50;
double SAFE_DISTANCE = 0.80;
double INHIBITED_RADIUS = 1.20;
double DIVERSITY_PATIENCE = 20;

std::string START_CONFIGURATION = "";
std::string GOAL_CONFIGURATION = "";
std::string SPACE_BOUNDS = "";
std::string MAP = "";
std::string OBJECT = "";
std::string GUIDING_OBJECT = "";
std::string CORRESPONDENCES = "";
std::string GUIDING_PATHS_DIRECTORY = "";
std::string GUIDING_PATHS_INDICES = "";
std::string OUTPUT = "";

bool CENTER = false;
bool SAVE_RESULTS = true;
/*	0 - generate guiding paths prior to the planning
	1 - load guiding paths prior to the planning
	2 - generate guiding paths only */
int MODE = 0;

// Transformation mapping the OBJECT to the GUIDING_OBJECT
double T[4][4];

namespace ob = ompl::base;
namespace og = ompl::geometric;

bool parse_args(int argc, char **argv);
const std::string exec(const char *cmd);

ob::PlannerPtr configuredRRTLIB(const ob::SpaceInformationPtr &si)
{
	og::RRTLIB *planner = new og::RRTLIB(si);

	planner->setRange(PLANNER_RANGE); // 0.0 results in automatic range evaluation
	planner->setGoalBias(GOAL_BIAS);
	planner->setPathBias(PATH_BIAS);
	planner->setSafeDistance(SAFE_DISTANCE);
	planner->setInhibitedRadius(INHIBITED_RADIUS);
	planner->setGuidingRadius(GUIDING_RADIUS);
	planner->setDiversityPatience(DIVERSITY_PATIENCE);

	planner->setGuidingPathsDirectory(GUIDING_PATHS_DIRECTORY);
	planner->setGuidingPathsIndices(GUIDING_PATHS_INDICES);

	planner->setGuidingPathsTransformation(T);

	planner->setMode(MODE);
	planner->setSaveResults(SAVE_RESULTS);
	planner->setOutput(OUTPUT);

	// Only used to save run info
	planner->setFilenames(MAP, GUIDING_OBJECT, OBJECT, CORRESPONDENCES);
	planner->setSetup(START_CONFIGURATION, GOAL_CONFIGURATION, SPACE_BOUNDS);

	return ompl::base::PlannerPtr(planner);
}

int main(int argc, char **argv)
{
	// Parse command line arguments
	if (!parse_args(argc, argv))
	{
		return EXIT_FAILURE;
	}

	// Prepare saving
	int runID = 0;
	if (SAVE_RESULTS && OUTPUT.empty())
	{
		// Determine save folder name and create it
		runID = std::stoi(exec("ls output/ | wc -l"));
		OUTPUT = "output/" + std::to_string(runID + 1) + "/";

		int res = std::system(("mkdir -p " + OUTPUT).c_str());
		if (res != 0)
		{
			std::cerr << "Couldn't create " << OUTPUT << "\n";
			exit(EXIT_FAILURE);
		}
	}

	// https://ompl.kavrakilab.org/geometricPlanningSE3.html
	// construct the state space we are planning in
	auto space(std::make_shared<ob::SE3StateSpace>());

	// Create a state space for the space we are planning in
	og::SimpleSetup ss(space);

	// Find transformation that maps the OBJECT onto the GUIDING_OBJECT
	// We can then transform the guiding paths computed for the GUIDING_OBJECT
	// to fit the OBJECT
	Matrix R(3, 3);
	Matrix t(3, 1);
	R.eye();
	t.zero();
	std::cout << "Calling ICP" << std::endl;
	getInitialTransformationICP(GUIDING_OBJECT, OBJECT, CORRESPONDENCES, R, t);

	std::cout << "Transformation results:"
			  << "\n";
	std::cout << "R:"
			  << "\n"
			  << R << "\n"
			  << "\n";
	std::cout << "t:"
			  << "\n"
			  << t << "\n"
			  << "\n";

	// Create a transformation matrix 'T'
	// from the rotation matrix 'R' and translation vector 't'
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			T[i][j] = R.val[i][j];
		}
	}
	for (int i = 0; i < 3; ++i)
	{
		T[i][3] = t.val[0][i];
	}
	T[3][0] = 0;
	T[3][1] = 0;
	T[3][2] = 0;
	T[3][3] = 1;

	const ob::StateValidityCheckerPtr vcp = std::make_shared<ValidityChecker>(ss.getSpaceInformation(),
																			  MAP, GUIDING_OBJECT, OBJECT);

	ss.setStateValidityChecker(vcp);
	ss.getSpaceInformation()->setStateValidityCheckingResolution(0.05);

	ob::RealVectorBounds bounds(3);
	const BoundingBox &bb = std::dynamic_pointer_cast<ValidityChecker>(vcp)->getObjectBoundingBox();

	// When the SPACE_BOUNDS string is not set, set the bounds to the
	// bounding box of the map, considering also the objects bounding box
	if (SPACE_BOUNDS.empty())
	{
		const BoundingBox &map_bb = std::dynamic_pointer_cast<ValidityChecker>(vcp)->getMapBoundingBox();

		bounds.setLow(0, map_bb.xmin - bb.xmin);
		bounds.setHigh(0, map_bb.xmax - bb.xmax);
		bounds.setLow(1, map_bb.ymin - bb.ymin);
		bounds.setHigh(1, map_bb.ymax - bb.ymax);
		bounds.setLow(2, map_bb.zmin - bb.zmin);
		bounds.setHigh(2, map_bb.zmax - bb.zmax);
	}
	else
	{
		std::stringstream space_bounds_str{SPACE_BOUNDS};
		double xmin, ymin, zmin, xmax, ymax, zmax;
		space_bounds_str >> xmin >> ymin >> zmin >> xmax >> ymax >> zmax;

		bounds.setLow(0, xmin);
		bounds.setHigh(0, xmax);
		bounds.setLow(1, ymin);
		bounds.setHigh(1, ymax);
		bounds.setLow(2, zmin);
		bounds.setHigh(2, zmax);
	}
	space->setBounds(bounds);

	// When CENTER is true, geometric center is calculated from the bounding box
	double x_c = 0, y_c = 0, z_c = 0;
	if (CENTER)
	{
		x_c = (bb.xmin + bb.xmax) / 2.0;
		y_c = (bb.ymin + bb.ymax) / 2.0;
		z_c = (bb.zmin + bb.zmax) / 2.0;
	}

	// Parse the START_CONFIGURATION and GOAL_CONFIGURATION strings,
	// which should be in the 'x y z rx ry rz' format
	std::stringstream start_str{START_CONFIGURATION};
	double start_x, start_y, start_z, start_qw, start_qx, start_qy, start_qz;
	start_str >> start_x >> start_y >> start_z >> start_qw >> start_qx >> start_qy >> start_qz;

	std::stringstream goal_str{GOAL_CONFIGURATION};
	double goal_x, goal_y, goal_z, goal_qw, goal_qx, goal_qy, goal_qz;
	goal_str >> goal_x >> goal_y >> goal_z >> goal_qw >> goal_qx >> goal_qy >> goal_qz;

	// Set START and GOAL configuration
	// Using StateSpacePtr patameter ('space') will automatically allocate
	// desired state. It still needs needs to be casted to allow setting the
	// values
	ob::ScopedState<> start(space);
	start->as<ob::SE3StateSpace::StateType>()->setXYZ(start_x - x_c, start_y - y_c, start_z - z_c);
	double start_d = std::sqrt(start_qw * start_qw + start_qx * start_qx + start_qy * start_qy + start_qz * start_qz);
	start->as<ob::SE3StateSpace::StateType>()->rotation().w = start_qw / start_d;
	start->as<ob::SE3StateSpace::StateType>()->rotation().x = start_qx / start_d;
	start->as<ob::SE3StateSpace::StateType>()->rotation().y = start_qy / start_d;
	start->as<ob::SE3StateSpace::StateType>()->rotation().z = start_qz / start_d;

	ob::ScopedState<> goal(space);
	goal->as<ob::SE3StateSpace::StateType>()->setXYZ(goal_x - x_c, goal_y - y_c, goal_z - z_c);
	double goal_d = std::sqrt(goal_qw * goal_qw + goal_qx * goal_qx + goal_qy * goal_qy + goal_qz * goal_qz);
	goal->as<ob::SE3StateSpace::StateType>()->rotation().w = goal_qw / goal_d;
	goal->as<ob::SE3StateSpace::StateType>()->rotation().x = goal_qx / goal_d;
	goal->as<ob::SE3StateSpace::StateType>()->rotation().y = goal_qy / goal_d;
	goal->as<ob::SE3StateSpace::StateType>()->rotation().z = goal_qz / goal_d;

	ss.setStartAndGoalStates(start, goal);

	ss.setPlanner(configuredRRTLIB(ss.getSpaceInformation()));
	ob::PlannerStatus solved = ss.solve(180.0);
	if (solved)
	{
		std::cout << "Found solution:" << std::endl;
		// print the path to screen
		ss.simplifySolution();
		ss.getSolutionPath().print(std::cout);
	}

	return 0;
}

bool parse_args(int argc, char **argv)
{
	InputParser input(argc, argv);

	if (input.cmdOptionExists("-h"))
	{
		std::cout << "--- RRT-LIB planner ---\n\n"
				  << "Usage: \n\n"
				  << "-M"
				  << "\t"
				  << "MODE"
				  << "\t\t\t"
				  << "Planner mode (default = 0)\n"
				  << "\t\t\t\t\t0 - generate guiding paths prior to the planning\n"
				  << "\t\t\t\t\t1 - load guiding paths prior to the planning\n"
				  << "\t\t\t\t\t2 - generate guiding paths only \n\n"
				  << "-s"
				  << "\t"
				  << "SAVE_RESULTS"
				  << "\t\t"
				  << "Set whether results should be saved (default = true)\n"
				  << "\n - FILES -\n"
				  << "-m"
				  << "\t"
				  << "MAP"
				  << "\t\t\t"
				  << "Map file (.off) location\n"
				  << "-o"
				  << "\t"
				  << "OBJECT"
				  << "\t\t\t"
				  << "Object file (.off) location\n"
				  << "-go"
				  << "\t"
				  << "GUIDING_OBJECT"
				  << "\t\t"
				  << "Guiding object file (.off) location\n"
				  << "-ct"
				  << "\t"
				  << "CENTER"
				  << "\t\t\t"
				  << "Whether the object should be centered on the start configuration (default = false)\n"
				  << "-c"
				  << "\t"
				  << "CORRESPONDENCES"
				  << "\t\t"
				  << "List of corresponding vertices between the object and the guiding object\n"
				  << "-pd"
				  << "\t"
				  << "PATHDIR"
				  << "\t\t\t"
				  << "Directory with guiding paths\n"
				  << "-p"
				  << "\t"
				  << "'P1,P2'"
				  << "\t\t\t"
				  << "Comma separated list of guiding path ids\n"
				  << "--out"
				  << "\t"
				  << "OUTPUT"
				  << "\t\t\t"
				  << "File to which the results will be saved\n"
				  << "\n - PROBLEM SETUP - \n"
				  << "-sc"
				  << "\t"
				  << "START_CONFIGURATION"
				  << "\t"
				  << "Start configuration 'x y z qw qx qy qz'\n"
				  << "-gc"
				  << "\t"
				  << "GOAL_CONFIGURATION"
				  << "\t"
				  << "Goal configuration 'x y z qw qx qy qz'\n"
				  << "-sb"
				  << "\t"
				  << "SPACE_BOUNDS"
				  << "\t\t"
				  << "The sampler bounds 'minx miny minz maxx maxy maxz'. Defaults to the map bounding box if not specified\n"
				  << "\n- RRT CONFIG -\n"
				  << "-gb"
				  << "\t"
				  << "GOAL_BIAS"
				  << "\t\t"
				  << "Goal bias (default = 0.05)\n"
				  << "-pb"
				  << "\t"
				  << "PATH_BIAS"
				  << "\t\t"
				  << "Probability of sampling along guiding paths (default = 0.80)\n"
				  << "-gr"
				  << "\t"
				  << "GUIDING_RADIUS"
				  << "\t\t"
				  << "Distance around the guiding path configurations where the samples are generated (default = 0.50)\n"
				  << "-sd"
				  << "\t"
				  << "SAFE_DISTANCE"
				  << "\t\t"
				  << "Distance around the start and goal configurations where configurations are not added to the inhibited regions (default = 0.80)\n"
				  << "-ir"
				  << "\t"
				  << "INHIBITED_RADIUS"
				  << "\t"
				  << "Distance around a configuration defining the inhibited region (default = 1.20)\n"
				  << "-dp"
				  << "\t"
				  << "DIVERSITY_PATIENCE"
				  << "\t"
				  << "How many paths need to be similar before stopping the search (default = 20)\n\n"
				  << "-h"
				  << "\t\t\t\t"
				  << "Show this help and exit"
				  << "\n";
		return false;
	}

	if (input.cmdOptionExists("-M"))
	{
		MODE = std::stoi(input.getCmdOption("-M"));
	}

	// FILES
	MAP = input.getCmdOption("-m");
	if (MAP.empty())
	{
		std::cout << "Path to a map file needs to be specified by the -m parameter\n";
		return false;
	}

	OBJECT = input.getCmdOption("-o");
	if (OBJECT.empty())
	{
		std::cout << "Path to an object file needs to be specified by the -o parameter\n";
		return false;
	}

	GUIDING_OBJECT = input.getCmdOption("-go");
	if (GUIDING_OBJECT.empty())
	{
		std::cout << "Path to the guiding object file needs to be specified by the -go parameter\n";
		return false;
	}

	CORRESPONDENCES = input.getCmdOption("-c");
	if (MODE != 2 && CORRESPONDENCES.empty())
	{
		std::cout << "Path to the file containing correspondences needs to be specified by the -c parameter\n";
		return false;
	}

	GUIDING_PATHS_DIRECTORY = input.getCmdOption("-pd");
	OUTPUT = input.getCmdOption("--out");

	if (input.cmdOptionExists("-ct"))
	{
		CENTER = std::stoi(input.getCmdOption("-ct"));
	}

	// Load all guiding paths if not specified otherwise
	GUIDING_PATHS_INDICES = input.getCmdOption("-p");
	if (MODE == 1 & GUIDING_PATHS_INDICES.empty())
	{
		for (const auto &entry : fs::directory_iterator(GUIDING_PATHS_DIRECTORY))
		{
			if (entry.path().extension() == ".path")
			{
				GUIDING_PATHS_INDICES += entry.path().stem().string() + ",";
			}
		}
		// Remove trailing comma
		GUIDING_PATHS_INDICES.pop_back();
		std::cout << "Guiding path indices not specified, loaded all: " << GUIDING_PATHS_INDICES << "\n";
	}

	// PROBLEM SETUP
	START_CONFIGURATION = input.getCmdOption("-sc");
	if (START_CONFIGURATION.empty())
	{
		std::cout << "Start configuration needs to be specified by the -sc parameter\n";
		return false;
	}

	GOAL_CONFIGURATION = input.getCmdOption("-gc");
	if (GOAL_CONFIGURATION.empty())
	{
		std::cout << "Goal configuration needs to be specified by the -gc parameter\n";
		return false;
	}

	SPACE_BOUNDS = input.getCmdOption("-sb");

	// RRT CONFIG
	if (input.cmdOptionExists("-gb"))
	{
		GOAL_BIAS = std::stod(input.getCmdOption("-gb"));
	}

	if (input.cmdOptionExists("-pr"))
	{
		PLANNER_RANGE = std::stod(input.getCmdOption("-pr"));
	}

	if (input.cmdOptionExists("-pb"))
	{
		PATH_BIAS = std::stoi(input.getCmdOption("-pb"));
	}

	if (input.cmdOptionExists("-gr"))
	{
		GUIDING_RADIUS = std::stod(input.getCmdOption("-gr"));
	}

	if (input.cmdOptionExists("-sd"))
	{
		SAFE_DISTANCE = std::stod(input.getCmdOption("-sd"));
	}

	if (input.cmdOptionExists("-ir"))
	{
		INHIBITED_RADIUS = std::stod(input.getCmdOption("-ir"));
	}

	if (input.cmdOptionExists("-dp"))
	{
		DIVERSITY_PATIENCE = std::stoi(input.getCmdOption("-dp"));
	}

	if (input.cmdOptionExists("-s"))
	{
		SAVE_RESULTS = std::stoi(input.getCmdOption("-s"));
	}

	return true;
}

// https://stackoverflow.com/questions/478898/how-do-i-execute-a-command-and-get-the-output-of-the-command-within-c-using-po
const std::string exec(const char *cmd)
{
	std::array<char, 128> buffer;
	std::string result;
	std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
	if (!pipe)
	{
		std::cerr << "popen() failed!";
		exit(EXIT_FAILURE);
	}
	while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
	{
		result += buffer.data();
	}
	return result;
}
