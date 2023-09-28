#ifndef __TYPES_EULER__
#define __TYPES_EULER__

#include <iostream>
#include <vector>
#include <string>
#include <limits>

#include "angle.h"
#include "nanoflann.hpp"

double get_random(const double from, const double to);
std::string align_to(std::string s, const unsigned int a, char c = ' ', bool ljust = false);
std::string align_to(double d, int a, bool ljust = false);

struct Point3
{
	Point3(const double ix = 0, const double iy = 0, const double iz = 0) : x(ix), y(iy), z(iz)
	{
	}
	Point3(const Point3 &obj) : x{obj.x}, y{obj.y}, z{obj.z}
	{
	}
	Point3 &operator=(const Point3 &rhs)
	{
		if (&rhs != this)
		{
			x = rhs.x;
			y = rhs.y;
			z = rhs.z;
		}
		return *this;
	}
	bool operator==(const Point3 &rhs) const
	{
		return (this->x == rhs.x) && (this->y == rhs.y) && (this->z == rhs.z);
	}
	bool operator!=(const Point3 &rhs) const
	{
		return !((*this) == rhs);
	}
	friend std::ostream &operator<<(std::ostream &os, const Point3 &r);
	double x, y, z;
};

struct Point3_3
{
	Point3_3(const double ix = 0, const double iy = 0, const double iz = 0, const double _alpha = 0, const double _beta = 0, const double _gamma = 0) : x(ix), y(iy), z(iz), alpha{_alpha}, beta{_beta}, gamma{_gamma}
	{
	}
	Point3_3(const Point3_3 &obj) : x{obj.x}, y{obj.y}, z{obj.z}, alpha{obj.alpha}, beta{obj.beta}, gamma{obj.gamma}
	{
	}
	Point3_3(const Point3 &obj) : x{obj.x}, y{obj.y}, z{obj.z}, alpha{0}, beta{0}, gamma{0}
	{
	}

	const std::vector<double> to_vec() const
	{
		return std::vector<double>{x, y, z, alpha, beta, gamma};
	}

	Point3_3 &operator=(const Point3_3 &rhs)
	{
		if (&rhs != this)
		{
			x = rhs.x;
			y = rhs.y;
			z = rhs.z;
			alpha = rhs.alpha;
			beta = rhs.beta;
			gamma = rhs.gamma;
		}
		return *this;
	}
	bool operator==(const Point3_3 &rhs) const
	{
		return (this->x == rhs.x) && (this->y == rhs.y) && (this->z == rhs.z) && (this->alpha == rhs.alpha) && (this->beta == rhs.beta) && (this->gamma == rhs.gamma);
	}
	bool operator!=(const Point3_3 &rhs) const
	{
		return !((*this) == rhs);
	}
	friend std::ostream &operator<<(std::ostream &os, const Point3_3 &r);

	double x, y, z;
	Angle alpha, beta, gamma;
};

struct Path
{
	Path() : path{}, id{0}, length{0}, time_to_find{0}
	{
	}
	std::vector<Point3_3> path;

	size_t id;
	size_t length;
	double time_to_find;
};

struct BoundingBox
{
	// For convenience, limits are initialized to inf values first
	// and calculated in loading functions
	BoundingBox()
	{
		xmin = std::numeric_limits<double>::max();
		ymin = std::numeric_limits<double>::max();
		zmin = std::numeric_limits<double>::max();

		xmax = std::numeric_limits<double>::min();
		ymax = std::numeric_limits<double>::min();
		zmax = std::numeric_limits<double>::min();
	}

	BoundingBox(const double _xmin, const double _xmax, const double _ymin, const double _ymax, const double _zmin, const double _zmax) : xmin{_xmin}, xmax{_xmax}, ymin{_ymin}, ymax{_ymax}, zmin{_zmin}, zmax{_zmax}
	{
	}

	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double zmin;
	double zmax;
};

// Holds RRT tree node, along with other info
struct Configuration
{
	Configuration(const Point3_3 &_point, const size_t _parent_index, const size_t _index = 0) : point{_point}, parent_index{_parent_index}, index{_index}
	{
	}
	const Point3_3 point;
	size_t parent_index;
	size_t index;
};

// Modified PointCloud from nanoflann/utils.h
struct KDTree
{
	std::vector<Configuration> pts;

	KDTree() = default;
	KDTree(const std::vector<Configuration> &_pts) : pts{_pts}
	{
	}
	// Must return the number of data points
	inline size_t kdtree_get_point_count() const
	{
		return pts.size();
	}

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	inline double kdtree_get_pt(const size_t idx, const size_t dim) const
	{
		if (dim == 0)
			return pts[idx].point.x;
		else if (dim == 1)
			return pts[idx].point.y;
		else if (dim == 2)
			return pts[idx].point.z;
		else if (dim == 3)
			return pts[idx].point.alpha;
		else if (dim == 4)
			return pts[idx].point.beta;
		else
			return pts[idx].point.gamma;
	}

	// For consistency and more readable code
	inline const Configuration &at(const size_t idx) const
	{
		return this->pts.at(idx);
	}

	inline size_t size() const
	{
		return this->pts.size();
	}

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX & /* bb */) const
	{
		return false;
	}
};

template <class T, class DataSource, typename _DistanceType = T>
struct My_Metric_Adaptor
{
	typedef T ElementType;
	typedef _DistanceType DistanceType;

	const DataSource &data_source;

	My_Metric_Adaptor(const DataSource &_data_source)
		: data_source(_data_source)
	{
	}

	inline DistanceType evalMetric(const T *a, const size_t b_idx,
								   size_t size) const
	{
		DistanceType result = DistanceType();
		for (size_t i = 0; i < size; ++i)
		{
			result += accum_dist(a[i], data_source.kdtree_get_pt(b_idx, i), i);
		}

		return result;
	}

	template <typename U, typename V>
	inline DistanceType accum_dist(const U a, const V b, const size_t size) const
	{
		if (size < 3)
		{
			return (a - b) * (a - b);
		}
		else
		{
			DistanceType diff = a - b;
			// From angle.h
			diff = wrap_to_pmPI(diff);
			return diff * diff;
		}
	}
};

typedef nanoflann::KDTreeSingleIndexDynamicAdaptor<
	My_Metric_Adaptor<double, KDTree>,
	KDTree,
	6 /* dim */
	>
	KDTree_index;

#endif // __TYPES_EULER__
