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

class polyhedron35 
{
	public:
	// initialize unit icosahedron
	polyhedron35();
	
	
	friend std::ostream& operator<< (std::ostream& os, const polyhedron35& p);

	void n_tesselate(size_t);
	void tesselate_to(size_t);

	private:
	void tesselate();
	std::vector<std::array<double,3>> vertices;
	std::vector<std::vector<size_t>> adj_mtx;
};

#endif
