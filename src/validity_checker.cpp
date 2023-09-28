#include <iostream>

#include <ompl/base/SpaceInformation.h>
#include <ompl/base/State.h> // https://ompl.kavrakilab.org/State_8h_source.html
#include <ompl/base/StateValidityChecker.h>
#include <ompl/base/spaces/SE3StateSpace.h> // https://ompl.kavrakilab.org/SE3StateSpace_8h_source.html
#include <ompl/base/spaces/SO3StateSpace.h> // https://ompl.kavrakilab.org/classompl_1_1base_1_1SO3StateSpace_1_1StateType.html
// https://github.com/ompl/ompl/blob/master/doc/markdown/optimalPlanningTutorial.md

#include "validity_checker.h"
#include "object.h"
#include "types.h"
#include "RAPID.H"

bool ValidityChecker::isValid(const ompl::base::State *state) const
{
    // Cast State to SE(3) representation (we know we work with this type
    // of space representation)
    const ompl::base::SE3StateSpace::StateType *state6D =
        state->as<ompl::base::SE3StateSpace::StateType>();

    // Get values:
    const ompl::base::SO3StateSpace::StateType &q = state6D->rotation();

    // Compute rotation matrix from quaternion
    double R[3][3];
    R[0][0] = 1 - 2 * q.y * q.y - 2 * q.z * q.z;
    R[0][1] = 2 * q.x * q.y - 2 * q.z * q.w;
    R[0][2] = 2 * q.x * q.z + 2 * q.y * q.w;

    R[1][0] = 2 * q.x * q.y + 2 * q.z * q.w;
    R[1][1] = 1 - 2 * q.x * q.x - 2 * q.z * q.z;
    R[1][2] = 2 * q.y * q.z - 2 * q.x * q.w;

    R[2][0] = 2 * q.x * q.z - 2 * q.y * q.w;
    R[2][1] = 2 * q.y * q.z + 2 * q.x * q.w;
    R[2][2] = 1 - 2 * q.x * q.x - 2 * q.y * q.y;

    // Get translation
    double T[3];
    T[0] = state6D->getX();
    T[1] = state6D->getY();
    T[2] = state6D->getZ();

    // Remove const by explicit casting to avoid compiler warnings
    RAPID_Collide(R, T, isSearchingPhase_ ? object_.model : guidingObject_.model,
                  (double(*)[3])map_.R, (double *)map_.T, map_.model,
                  RAPID_FIRST_CONTACT);

    return RAPID_num_contacts == 0;
}

const BoundingBox &ValidityChecker::getObjectBoundingBox() const
{
    return object_.bb;
}

const BoundingBox &ValidityChecker::getMapBoundingBox() const
{
    return map_.bb;
}
