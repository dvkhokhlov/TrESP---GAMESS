#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

#include <iostream>
#include <iomanip>
#include <cassert>
#include <string>
#include <stdexcept>
#include <vector>
#include <memory>
#include <limits>
#include <ctgmath>
#include <array>
#include <algorithm>

#include "algebra.h"

struct edge
{
	bool dividedQ;
	size_t spawned_v;
	std::array<size_t,2> v;
};

struct face
{	
	std::array<size_t, 3> e;
};

class polyhedron35 
{
	public:
	
	// initialize unit icosahedron
	polyhedron35();
	
	friend std::ostream& operator<< (std::ostream& os, const polyhedron35& p);
	
	void tesselate();
	
	private:
	inline double powm1 (double a, int pow)
	{
		return pow%2 ? -a : a;
	}
	
	void gen_edges();
	void gen_faces();
	
	const size_t nvert0 = 12;
	const size_t nedge0 = 30;
	const size_t nface0 = 20;
	const double unit_len2 = 4.;
	const double t = (1. + sqrt(5.))/2.;
	const double r = sqrt(1. + t*t);
	
	std::vector<std::array<double,3>> vertices;
	std::vector<edge> edges;
	std::vector<face> faces;
};

#endif
