#include "ompl/base/spaces/SE3StateSpace.h"
#include "ompl/base/spaces/SO3StateSpace.h"

#include <cmath>

namespace ob = ompl::base;

void TtoSE3(const double T[4][4], ob::SE3StateSpace::StateType *state)
{
	// https://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
	double trace = T[0][0] + T[1][1] + T[2][2];
	if (trace > 0)
	{
		double s = 0.5f / std::sqrt(trace + 1.0f);
		state->rotation().w = 0.25f / s;
		state->rotation().x = (T[2][1] - T[1][2]) * s;
		state->rotation().y = (T[0][2] - T[2][0]) * s;
		state->rotation().z = (T[1][0] - T[0][1]) * s;
	}
	else
	{
		if (T[0][0] > T[1][1] && T[0][0] > T[2][2])
		{
			double s = 2.0f * std::sqrt(1.0f + T[0][0] - T[1][1] - T[2][2]);
			state->rotation().w = (T[2][1] - T[1][2]) / s;
			state->rotation().x = 0.25f * s;
			state->rotation().y = (T[0][1] + T[1][0]) / s;
			state->rotation().z = (T[0][2] + T[2][0]) / s;
		}
		else if (T[1][1] > T[2][2])
		{
			double s = 2.0f * std::sqrt(1.0f + T[1][1] - T[0][0] - T[2][2]);
			state->rotation().w = (T[0][2] - T[2][0]) / s;
			state->rotation().x = (T[0][1] + T[1][0]) / s;
			state->rotation().y = 0.25f * s;
			state->rotation().z = (T[1][2] + T[2][1]) / s;
		}
		else
		{
			double s = 2.0f * std::sqrt(1.0f + T[2][2] - T[0][0] - T[1][1]);
			state->rotation().w = (T[1][0] - T[0][1]) / s;
			state->rotation().x = (T[0][2] + T[2][0]) / s;
			state->rotation().y = (T[1][2] + T[2][1]) / s;
			state->rotation().z = 0.25f * s;
		}
	}

	state->setXYZ(T[0][3], T[1][3], T[2][3]);
}

void M4x4(const double A[4][4], const double B[4][4], double C[4][4])
{
	// C = A * B
	double num;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			num = 0;
			for (int k = 0; k < 4; k++)
			{
				num += A[i][k] * B[k][j];
			}
			C[i][j] = num;
		}
	}
}

void quatMultiply(const ob::SO3StateSpace::StateType *q1, const ob::SO3StateSpace::StateType *q2, ob::SO3StateSpace::StateType *q)
{
	double x, y, z, w;
	x = q1->x * q2->w + q1->y * q2->z - q1->z * q2->y + q1->w * q2->x;
	y = -q1->x * q2->z + q1->y * q2->w + q1->z * q2->x + q1->w * q2->y;
	z = q1->x * q2->y - q1->y * q2->x + q1->z * q2->w + q1->w * q2->z;
	w = -q1->x * q2->x - q1->y * q2->y - q1->z * q2->z + q1->w * q2->w;

	q->x = x;
	q->y = y;
	q->z = z;
	q->w = w;
}
