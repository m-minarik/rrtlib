#ifndef __VALIDITY_CHECKER_H__
#define __VALIDITY_CHECKER_H__

#include <ompl/base/SpaceInformation.h>
#include <ompl/base/State.h> // https://ompl.kavrakilab.org/State_8h_source.html
#include <ompl/base/StateValidityChecker.h>

#include "object.h"
#include "icpPointToPlane.h"

class ValidityChecker : public ompl::base::StateValidityChecker
{
public:
	ValidityChecker(const ompl::base::SpaceInformationPtr &si, const std::string &map_filename,
					const std::string &guiding_object_filename, const std::string &object_filename,
					const Matrix *R = nullptr, const Matrix *t = nullptr) : ompl::base::StateValidityChecker(si), map_(map_filename), guidingObject_(guiding_object_filename), object_(object_filename, R, t)

	/* Initializes validity checker made for RRT-LIB
	 *
	 * R, t - initial transformation of the object (mapping object_ to guidingObject_)
	 */
	{
	}

	virtual bool isValid(const ompl::base::State *state) const;

	void startGuidingPhase()
	{
		isSearchingPhase_ = false;
	}
	void startSearchPhase()
	{
		isSearchingPhase_ = true;
	}

	const BoundingBox &getObjectBoundingBox() const;
	const BoundingBox &getMapBoundingBox() const;

private:
	Object map_;

	Object guidingObject_;
	Object object_;

	bool isSearchingPhase_{true};
};

#endif //__VALIDITY_CHECKER_H__
