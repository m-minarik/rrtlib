#ifndef __OMPL_UTILS__
#define __OMPL_UTILS__

#include "matrix.h" // from libcip
#include "ompl/base/spaces/SE3StateSpace.h"
#include "ompl/base/spaces/SO3StateSpace.h"

#include <cmath>

namespace ob = ompl::base;

void TtoSE3(const double T[4][4], ob::SE3StateSpace::StateType *state);
void quatMultiply(const ob::SO3StateSpace::StateType *q1, const ob::SO3StateSpace::StateType *q2, ob::SO3StateSpace::StateType *q);
void M4x4(const double A[4][4], const double B[4][4], double C[4][4]);

#endif
