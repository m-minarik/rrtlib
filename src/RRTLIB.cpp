/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2023, Michal Minařík
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *********************************************************************/

#include <limits>

#include <ompl/geometric/PathGeometric.h>
#include <ompl/base/goals/GoalSampleableRegion.h>
#include <ompl/geometric/planners/rrt/RRTLIB.h>
#include <ompl/tools/config/SelfConfig.h>
#include <ompl/base/goals/GoalState.h>

#include <ompl/base/SpaceInformation.h>
#include <ompl/base/State.h>				// https://ompl.kavrakilab.org/State_8h_source.html
#include <ompl/base/spaces/SE3StateSpace.h> // https://ompl.kavrakilab.org/SE3StateSpace_8h_source.html
#include <ompl/base/spaces/SO3StateSpace.h> // https://ompl.kavrakilab.org/classompl_1_1base_1_1SO3StateSpace_1_1StateType.html
// https://github.com/ompl/ompl/blob/master/doc/markdown/optimalPlanningTutorial.md

#include "validity_checker.h"
#include "ompl_utils.h"

// Saving
#include <iostream>
#include <iomanip>
#include <fstream>
#include <filesystem>

#include <chrono>
typedef std::chrono::time_point<std::chrono::steady_clock> time_point;

#include "json.hpp" // https://github.com/nlohmann/json.git
namespace fs = std::filesystem;

void ompl::geometric::RRTLIB::saveResults(const std::vector<ompl::geometric::PathGeometricPtr> &paths) const
{
	// Convert relative paths to absolute to make rendering easier
	nlohmann::json j = {
		{"params", {{"map", fs::absolute(fs::path{fileMap_})},

					{"object", fileObject_.empty() ? "" : fs::absolute(fs::path{fileObject_})},
					{"guiding_object", fileGuidingObject_.empty() ? "" : fs::absolute(fs::path{fileGuidingObject_})},
					{"correspondences", fileCorrespondences_.empty() ? "" : fs::absolute(fs::path{fileCorrespondences_})},
					{"start", startConfigurationStr_},
					{"goal", goalConfigurationStr_},
					{"space_bounds", spaceBoundsStr_},
					{"range", maxDistance_},
					{"goal_bias", goalBias_},
					{"path_bias", pathBias_},
					{"safe_distance", safeDistance_},
					{"inhibited_radius", inhibitedRadius_},
					{"diversity_patience", diversityPatience_},
					{"guiding_paths_indices", guidingPathsIndices_},
					{"guiding_paths_directory", guidingPathsDirectory_}}},
		{"paths", nullptr}};

	// Saves each path in format
	// number of steps
	// x y z qx qy qz qw
	int path_id = 0;
	for (const auto &path : paths)
	{
		std::string filename = output_ + "/" + std::to_string(path_id) + ".path";
		std::ofstream output(filename, std::ofstream::out);
		output << std::fixed << std::setprecision(3);

		output << path->getStateCount() << "\n";
		path->printAsMatrix(output);
		output.close();

		j["paths"][std::to_string(path_id)] = {
			{"length", path->getStateCount()},
			{"path_filename", fs::absolute(fs::path{filename})}};

		++path_id;
	}

	// write prettified JSON
	std::ofstream o(output_ + "/info.json");
	o << std::setw(4) << j << std::endl;
	o.close();

	// Create a symlink to the object file
	system(("ln -s " + fs::absolute(fs::path{fileObject_}).string() + " " + output_ + "/object.off").c_str());
}

#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif

ompl::geometric::RRTLIB::RRTLIB(const base::SpaceInformationPtr &si, bool addIntermediateStates)
	: base::Planner(si, addIntermediateStates ? "RRTLIBintermediate" : "RRTLIB")
{
	specs_.approximateSolutions = true;
	specs_.directed = true;

	Planner::declareParam<double>("range", this, &RRTLIB::setRange, &RRTLIB::getRange, "0.:1.:10000.");
	Planner::declareParam<double>("goal_bias", this, &RRTLIB::setGoalBias, &RRTLIB::getGoalBias, "0.:.05:1.");
	Planner::declareParam<double>("path_bias", this, &RRTLIB::setPathBias, &RRTLIB::getPathBias, "0.:.01:1.");
	Planner::declareParam<double>("safe_distance", this, &RRTLIB::setSafeDistance, &RRTLIB::getSafeDistance, "0.:.05:10.");
	Planner::declareParam<double>("inhibited_radius", this, &RRTLIB::setInhibitedRadius, &RRTLIB::getInhibitedRadius, "0.:.05:10.");
	Planner::declareParam<double>("guiding_radius", this, &RRTLIB::setGuidingRadius, &RRTLIB::getGuidingRadius, "0.:.05:10.");
	Planner::declareParam<int>("diversity_patience", this, &RRTLIB::setDiversityPatience, &RRTLIB::getDiversityPatience, "0:1:100");
	Planner::declareParam<bool>("intermediate_states", this, &RRTLIB::setIntermediateStates, &RRTLIB::getIntermediateStates, "0,1");
	Planner::declareParam<std::string>("guiding_paths_directory", this, &RRTLIB::setGuidingPathsDirectory, &RRTLIB::getGuidingPathsDirectory);
	Planner::declareParam<std::string>("guiding_paths_indices", this, &RRTLIB::setGuidingPathsIndices, &RRTLIB::getGuidingPathsIndices);
	Planner::declareParam<bool>("save_results", this, &RRTLIB::setSaveResults, &RRTLIB::getSaveResults, "0,1");
	Planner::declareParam<int>("mode", this, &RRTLIB::setMode, &RRTLIB::getMode, "0,1,2");

	addIntermediateStates_ = addIntermediateStates;
}

ompl::geometric::RRTLIB::~RRTLIB()
{
	freeMemory();
}

void ompl::geometric::RRTLIB::clear()
{
	Planner::clear();
	sampler_.reset();
	freeMemory();

	guidingPaths_.clear();
	lastGoalMotion_ = nullptr;
}

void ompl::geometric::RRTLIB::setup()
{
	Planner::setup();
	tools::SelfConfig sc(si_, getName());
	sc.configurePlannerRange(maxDistance_);

	if (!nn_)
		nn_.reset(tools::SelfConfig::getDefaultNearestNeighbors<Motion *>(this));
	nn_->setDistanceFunction([this](const Motion *a, const Motion *b)
							 { return distanceFunction(a, b); });

	if (!inhibitedNN_)
		inhibitedNN_.reset(tools::SelfConfig::getDefaultNearestNeighbors<InhibitedRegion *>(this));
	inhibitedNN_->setDistanceFunction([this](const InhibitedRegion *a, const InhibitedRegion *b)
									  { return distanceFunction(a, b); });
}

void ompl::geometric::RRTLIB::freeMemory()
{
	freeNN_();
	freeInhibitedNN_();
}

void ompl::geometric::RRTLIB::freeNN_()
{
	if (nn_)
	{
		std::vector<Motion *> motions;
		nn_->list(motions);
		for (auto &motion : motions)
		{
			if (motion->state != nullptr)
				si_->freeState(motion->state);
			delete motion;
		}
		nn_->clear();
	}
}

void ompl::geometric::RRTLIB::freeInhibitedNN_()
{
	if (inhibitedNN_)
	{
		std::vector<InhibitedRegion *> regions;
		inhibitedNN_->list(regions);
		for (auto &region : regions)
		{
			if (region->state != nullptr)
				si_->freeState(region->state);
			delete region;
		}
		inhibitedNN_->clear();
	}
}

ompl::base::PlannerStatus ompl::geometric::RRTLIB::solve(const base::PlannerTerminationCondition &ptc)
{
	// There will be a single start, so get 0th element
	const ompl::base::State *start = pdef_->getStartState(0);
	const ompl::base::State *goal = pdef_->getGoal()->as<ompl::base::GoalState>()->getState();

	// Load specified guiding paths
	if (mode_ == 1)
	{
		if (guidingPathsDirectory_.empty())
		{
			std::cerr << "Guiding path directory is not set!\n";
			exit(EXIT_FAILURE);
		}
		std::stringstream pathIDs{guidingPathsIndices_};
		std::string pathID;
		size_t i = 0;
		while (std::getline(pathIDs, pathID, ','))
		{
			std::string pathFilename = guidingPathsDirectory_ + pathID + ".path";
			// From RRT_LIB
			if (loadGuidingPath(pathFilename, guidingPaths_))
			{
				for (size_t j = 0; j < guidingPaths_.back()->getStateCount(); ++j)
				{
					const ompl::base::State *s = guidingPaths_.back()->getState(j);

					double startDist = si_->distance(s, start);
					double goalDist = si_->distance(s, goal);

					if (startDist > safeDistance_ && goalDist > safeDistance_)
					{
						auto *ir = new InhibitedRegion(si_, i, j);
						si_->copyState(ir->state, s);
						inhibitedNN_->add(ir);
					}
				}
				i += 1;
			}
		}
	}
	else
	{
		// Generate guiding paths
		int n_same = 0;
		size_t path_index = 0;

		// Uses guidingObject_ to search for guiding paths
		dynamic_cast<ValidityChecker &>(*si_->getStateValidityChecker()).startGuidingPhase();

		while (n_same < diversityPatience_)
		{
			ompl::base::PlannerStatus status = findPath(ptc, {});
			if (status)
			{
				double path_distinctness = pathDistinctness(lastFoundPath_);
				OMPL_INFORM("%s: Found next path with distinctness %f", getName().c_str(), path_distinctness);

				if (path_distinctness > inhibitedRadius_)
				{
					OMPL_INFORM("%s: Path added to guiding paths: %lu", getName().c_str(), guidingPaths_.size());
					guidingPaths_.push_back(lastFoundPath_);
					n_same = 0;
				}
				else
				{
					OMPL_INFORM("%s: Path not added to guiding paths", getName().c_str());
					n_same += 1;
				}

				// If path was found, update inhibited regions
				for (size_t j = 0; j < lastFoundPath_->getStateCount(); ++j)
				{
					const ompl::base::State *s = lastFoundPath_->getState(j);

					double startDist = si_->distance(s, start);
					double goalDist = si_->distance(s, goal);

					if (startDist > safeDistance_ && goalDist > safeDistance_)
					{
						auto *ir = new InhibitedRegion(si_, path_index, j);
						si_->copyState(ir->state, s);
						inhibitedNN_->add(ir);
					}
				}
				path_index += 1;
			}
			else
			{
				// If findPath ended due to termination condition, break
				break;
			}

			/* Clear only planning data (nn_) not inhibited regions (inhibitedNN_ and inhibitedRegions_) */
			Planner::clear();
			freeNN_();
		}

		// mode == 2: generate guiding paths only
		if (saveResults_ && mode_ == 2)
		{
			saveResults(guidingPaths_);
			pdef_->addSolutionPath(lastFoundPath_, false, 0.0, getName());
			return {true, false};
		}
	}

	Planner::clear();
	freeNN_();
	freeInhibitedNN_();

	// Uses original object to search for a path
	dynamic_cast<ValidityChecker &>(*si_->getStateValidityChecker()).startSearchPhase();

	ompl::base::PlannerStatus status = findPath(ptc, guidingPaths_);
	if (status)
	{
		// `approximate` would be false and `approxdif` 0.0 anyways in the case of
		// addSolutionPath being called, so there is no need to somehow return it
		// from findPath
		// pdef_->addSolutionPath(lastFoundPath_, approximate, approxdif, getName());
		OMPL_INFORM("%s: Found final path", getName().c_str());
		pdef_->addSolutionPath(lastFoundPath_, false, 0.0, getName());

		if (saveResults_)
		{
			saveResults({lastFoundPath_});
		}
	}
	return status;
}

double ompl::geometric::RRTLIB::getSimilarity(const std::shared_ptr<ompl::geometric::PathGeometric> &a, const std::shared_ptr<ompl::geometric::PathGeometric> &b) const
{
	double similarity = 0;

	for (size_t i = 0; i < a->getStateCount(); ++i)
	{
		// For each point in the first path, find nearest point in the second
		// path and add distance between them to the `similarity`
		similarity += si_->distance(a->getState(i), b->getState(b->getClosestIndex(a->getState(i))));
	}
	return similarity / a->getStateCount();
}

double ompl::geometric::RRTLIB::pathDistinctness(const std::shared_ptr<ompl::geometric::PathGeometric> &new_path) const
{
	double distinctness = std::numeric_limits<double>::max();
	for (size_t i = 0; i < guidingPaths_.size(); ++i)
	{
		double s_1 = getSimilarity(guidingPaths_[i], new_path);
		double s_2 = getSimilarity(new_path, guidingPaths_[i]);
		double s = MAX(s_1, s_2);
		distinctness = MIN(distinctness, s);
	}

	return distinctness;
}

ompl::base::PlannerStatus ompl::geometric::RRTLIB::findPath(const base::PlannerTerminationCondition &ptc, const std::vector<ompl::geometric::PathGeometricPtr> &guidingPaths)
{
	checkValidity();
	base::Goal *goal = pdef_->getGoal().get();
	auto *goal_s = dynamic_cast<base::GoalSampleableRegion *>(goal);

	while (const base::State *st = pis_.nextStart())
	{
		auto *motion = new Motion(si_);
		si_->copyState(motion->state, st);
		nn_->add(motion);
	}

	if (nn_->size() == 0)
	{
		OMPL_ERROR("%s: There are no valid initial states!", getName().c_str());
		return base::PlannerStatus::INVALID_START;
	}

	if (!sampler_)
		sampler_ = si_->allocStateSampler();

	OMPL_INFORM("%s: Starting planning with %u states already in datastructure", getName().c_str(), nn_->size());

	// Compute sizes beforehand
	size_t G_size = guidingPaths.size();
	size_t R_size = inhibitedNN_->size();

	// Initialize RRT_LIB variables
	std::vector<size_t> active_waypoints(G_size, 0);

	// Prepare inhibited regions
	std::vector<std::vector<size_t>> attempts{};
	std::vector<size_t> a_num{};
	std::vector<InhibitedRegion *> regions;
	inhibitedNN_->list(regions);
	for (const auto &region : regions)
	{
		while (attempts.size() < region->pathIndex + 1)
		{
			attempts.emplace_back();
			a_num.push_back(0);
		}
		while (attempts[region->pathIndex].size() < region->regionIndex)
		{
			attempts[region->pathIndex].push_back(0);
		}
	}

	size_t attempts_sum = 0;

	// Create variables beforehand to increase execution speed
	Point3_3 sampled_point, q_ij;
	size_t nearest_index;

	double dx, dy, dz, dalpha, dbeta, dgamma;
	unsigned int step, substep, number_of_steps;
	Point3_3 new_point, temp_point;
	bool collision;
	size_t i, j, k;
	size_t parent_index;

	double p;

	/* inhibitedNN_->nearest() needs inhibitedRegion instance, create it and
	   allocate state once */
	auto *tempregion = new InhibitedRegion(si_, 0, 0);
	InhibitedRegion *nregion = nullptr;

	Motion *solution = nullptr;
	Motion *approxsol = nullptr;
	double approxdif = std::numeric_limits<double>::infinity();
	auto *rmotion = new Motion(si_);
	base::State *rstate = rmotion->state;
	base::State *xstate = si_->allocState();

	std::pair<ompl::base::State *, double> lastValid;
	lastValid.first = si_->allocState();

	while (!ptc)
	{
		/* sample random state (with waypoint and goal biasing) */
		if ((goal_s != nullptr) && rng_.uniform01() < goalBias_ && goal_s->canSample())
		{
			// Sample around goal state
			goal_s->sampleGoal(rstate);
		}
		else if (G_size && rng_.uniform01() < pathBias_)
		{
			// Get random guiding path index
			i = rng_.uniformInt(0, G_size - 1);

			// Sample a state near the current waypoint at path i
			sampler_->sampleGaussian(rstate, guidingPaths[i]->getState(active_waypoints[i]), guidingRadius_);
		}
		else
		{
			// Sample random state inside the sampling space
			sampler_->sampleUniform(rstate);
		}

		/* find closest state in the tree */
		Motion *nmotion = nn_->nearest(rmotion);

		base::State *dstate = rstate;

		/* find state to add */
		/* If distance is larger than maxDistance_, interpolate new state
		   between nearest state and sampled state */
		double d = si_->distance(nmotion->state, rstate);
		if (d > maxDistance_)
		{
			si_->getStateSpace()->interpolate(nmotion->state, rstate, maxDistance_ / d, xstate);
			dstate = xstate;
		}

		/* Calls ValidityChecker.isValid(state) for states between
		   nmotion->state and dstate in this case */
		if (si_->checkMotion(nmotion->state, dstate))
		{
			p = 1;
			if (R_size)
			{
				/* Find nearest inhibited region state */
				si_->copyState(tempregion->state, nmotion->state);
				nregion = inhibitedNN_->nearest(tempregion);

				i = nregion->pathIndex;
				j = nregion->regionIndex;

				if (si_->distance(nregion->state, tempregion->state) < inhibitedRadius_)
				{
					attempts[i][j] += 1;
					attempts_sum += 1;
					p = calculateP_(i, j, attempts, attempts_sum);
					}
			}

			if (rng_.uniform01() < p)
			{
				/* Add motion to the tree */
				auto *motion = new Motion(si_);
				si_->copyState(motion->state, dstate);
				motion->parent = nmotion;
				nn_->add(motion);

				nmotion = motion;

				/* Check active waypoint on each path */
				for (i = 0; i < G_size; ++i)
				{
					while (si_->distance(guidingPaths[i]->getState(active_waypoints[i]), nmotion->state) < guidingRadius_)
					{
						if (active_waypoints[i] < guidingPaths[i]->getStateCount() - 1)
						{
							active_waypoints[i] += 1;
						}
						else
						{
							break;
						}
					}
				}
			}
			else
			{
				continue;
			}
			double dist = 0.0;
			bool sat = goal->isSatisfied(nmotion->state, &dist);
			if (sat)
			{
				approxdif = dist;
				solution = nmotion;
				break;
			}
			if (dist < approxdif)
			{
				approxdif = dist;
				approxsol = nmotion;
			}
		}
	}

	bool solved = false;
	bool approximate = false;
	if (solution == nullptr)
	{
		solution = approxsol;
		approximate = true;
	}
	else
	{
		solved = true;
	}

	// build the solution path
	if (solution != nullptr)
	{
		lastGoalMotion_ = solution;

		// construct the solution path
		std::vector<Motion *> mpath;
		while (solution != nullptr)
		{
			mpath.push_back(solution);
			solution = solution->parent;
		}

		// set the solution path
		lastFoundPath_ = std::make_shared<PathGeometric>(si_);
		for (int i = mpath.size() - 1; i >= 0; --i)
			lastFoundPath_->append(mpath[i]->state);
	}

	si_->freeState(xstate);

	if (rmotion->state != nullptr)
		si_->freeState(rmotion->state);
	delete rmotion;

	if (tempregion->state != nullptr)
		si_->freeState(tempregion->state);
	delete tempregion;

	if (lastValid.first != nullptr)
		si_->freeState(lastValid.first);

	return {solved, approximate};
}

void ompl::geometric::RRTLIB::getPlannerData(base::PlannerData &data) const
{
	Planner::getPlannerData(data);

	std::vector<Motion *> motions;
	if (nn_)
		nn_->list(motions);

	if (lastGoalMotion_ != nullptr)
		data.addGoalVertex(base::PlannerDataVertex(lastGoalMotion_->state));

	for (auto &motion : motions)
	{
		if (motion->parent == nullptr)
			data.addStartVertex(base::PlannerDataVertex(motion->state));
		else
			data.addEdge(base::PlannerDataVertex(motion->parent->state), base::PlannerDataVertex(motion->state));
	}
}

inline double ompl::geometric::RRTLIB::calculateP_(const size_t &i, const size_t &j, const std::vector<std::vector<size_t>> &attempts, const size_t &attempts_sum) const
{
	size_t A_ij = 0;
	size_t B_ij = 0;
	for (size_t k = 0; k < attempts[i].size(); ++k)
	{
		if (k > j)
		{
			A_ij = MAX(A_ij, attempts[i][k]);
		}
		else
		{
			B_ij = MAX(B_ij, attempts[i][k]);
		}
	}

	if (A_ij == 0)
	{
		return exp(-(double)B_ij / (double)attempts_sum);
	}
	else
	{
		return 0.0;
	}
}

bool ompl::geometric::RRTLIB::loadGuidingPath(const std::string &pathFilename, std::vector<ompl::geometric::PathGeometricPtr> &guidingPaths)
{
	std::ifstream input(pathFilename);
	if (!input.is_open())
	{
		std::cout << "Cannot open " << pathFilename << "\n";
		return false;
	}

	// Coordinates
	double x, y, z, qx, qy, qz, qw;
	size_t pathLength;

	std::string line;
	size_t lineNo = 1;

	std::shared_ptr<PathGeometric> path = std::make_shared<PathGeometric>(si_);

	// Transformation matrix
	double T_orig[4][4];
	double T_new[4][4];

	while (std::getline(input, line))
	{
		std::stringstream stream(line);
		if (lineNo == 1)
		{
			stream >> pathLength;
		}
		else if (lineNo - 2 < pathLength)
		{
			base::State *state = si_->allocState();
			stream >> x >> y >> z >> qx >> qy >> qz >> qw;

			// Compute rotation matrix from quaternion
			T_orig[0][0] = 1 - 2 * qy * qy - 2 * qz * qz;
			T_orig[0][1] = 2 * qx * qy - 2 * qz * qw;
			T_orig[0][2] = 2 * qx * qz + 2 * qy * qw;

			T_orig[1][0] = 2 * qx * qy + 2 * qz * qw;
			T_orig[1][1] = 1 - 2 * qx * qx - 2 * qz * qz;
			T_orig[1][2] = 2 * qy * qz - 2 * qx * qw;

			T_orig[2][0] = 2 * qx * qz - 2 * qy * qw;
			T_orig[2][1] = 2 * qy * qz + 2 * qx * qw;
			T_orig[2][2] = 1 - 2 * qx * qx - 2 * qy * qy;

			// Get translation
			T_orig[0][3] = x;
			T_orig[1][3] = y;
			T_orig[2][3] = z;
			T_orig[3][0] = 0;
			T_orig[3][1] = 0;
			T_orig[3][2] = 0;
			T_orig[3][3] = 1;

			// T_new = T_orig * T_
			M4x4(T_orig, T_, T_new);

			ompl::base::SE3StateSpace::StateType *stateSE3 = state->as<ompl::base::SE3StateSpace::StateType>();

			// Convert transformation matrix T to a SE3 state
			TtoSE3(T_new, stateSE3);

			path->append(state);
		}
		++lineNo;
	}
	input.close();

	guidingPaths.push_back(path);
	return true;
}
