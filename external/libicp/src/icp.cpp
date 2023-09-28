/*
   Copyright 2011. All rights reserved.
   Institute of Measurement and Control Systems
   Karlsruhe Institute of Technology, Germany

   Authors: Andreas Geiger

   libicp is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation; either version 2 of the License, or any later version.

   libicp is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
   PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with
   libicp; if not, write to the Free Software Foundation, Inc., 51 Franklin
   Street, Fifth Floor, Boston, MA 02110-1301, USA
 */

#include "icp.h"
#include <chrono>
#include <fstream>

typedef std::chrono::time_point<std::chrono::steady_clock> time_point;

using namespace std;

Icp::Icp (double *M, const int32_t M_num, const int32_t dim, const std::string &save_dir) :
	dim(dim), save_dir(save_dir), max_iter(15), min_delta(1e-10)
{

	  // check for correct dimensionality
	if (dim != 2 && dim != 3) {
		cout << "ERROR: LIBICP works only for data of dimensionality 2 or 3" << endl;
		M_tree = 0;
		return;
	}

	  // check for minimum number of points
	if (M_num < 5) {
		cout << "ERROR: LIBICP works only with at least 5 model points" << endl;
		M_tree = 0;
		return;
	}

	  // copy model points to M_data
	M_data.resize(boost::extents[M_num][dim]);
	for (int32_t m = 0; m < M_num; m++)
		for (int32_t n = 0; n < dim; n++)
			M_data[m][n] = (float)M[m * dim + n];

	  // build a kd tree from the model point cloud
	M_tree = new kdtree::KDTree(M_data);
}

Icp::~Icp ()
{
	if (M_tree)
		delete M_tree;
}

double Icp::fit (double *T, const int32_t T_num, Matrix &R, Matrix &t, const double indist)
{

	  // make sure we have a model tree
	if (!M_tree) {
		cout << "ERROR: No model available." << endl;
		return -1;
	}

	  // check for minimum number of points
	if (T_num < 5) {
		cout << "ERROR: Icp works only with at least 5 template points" << endl;
		return -1;
	}

	  // set active points
	vector<int32_t> active;
	if (indist <= 0) {
		active.clear();
		for (int32_t i = 0; i < T_num; i++)
			active.push_back(i);
	} else {
		active = getInliers(T, T_num, R, t, indist);
	}

	  // run icp
	return fitIterate(T, T_num, R, t, active);
}

double Icp::fitIterate(double *T, const int32_t T_num, Matrix &R, Matrix &t, const std::vector<int32_t> &active)
{

	  // check if we have at least 5 active points
	if (active.size() < 5)
		return -1;

	// Time measurement
	// time_point start_t = std::chrono::steady_clock::now();
	// std::ofstream time_stats(save_dir + "time.csv", std::ios::app);
	// std::ofstream metric_stats(save_dir + "metric.csv", std::ios::app);

	  // iterate until convergence
	double d;
	for (int32_t iter = 0; iter < max_iter; iter++) {
		d = fitStep(T, T_num, R, t, active);

		// Get milliseconds resolution
		// time_point end_t = std::chrono::steady_clock::now();
		// double ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - start_t).count();

		  // time_stats << iter + 1 << "," << ms / 1000.0 << "\n";
		  // metric_stats << iter + 1 << "," << d << "\n";
		  // std::cout << iter + 1 << "," << d << "\n";
		if (d < min_delta)
			break;
	}

	return d;
	// time_stats.close();
	// metric_stats.close();
}
