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
#include "chemistry.h"

class polyhedron35 
{
	public:
	polyhedron35();
	polyhedron35(size_t);
	
	friend std::ostream& operator<< (std::ostream& os, const polyhedron35& p);

	void n_tesselate(size_t);
	void tesselate_to(size_t);
	const std::vector<std::array<double,3>>& get_vertices(); 

	private:
	void tesselate();
	std::vector<std::array<double,3>> vertices;
	std::vector<std::vector<size_t>> adj_mtx;
	
	public:
	const double t = (1. + sqrt(5.))/2.;
	const double unit_len2 = 4.;
	const double r_spher = sqrt(1. + t*t);
};

class grid
{
	public:
	grid(const std::vector<std::array<double,3>>& sphere, double r_spher);
	grid(const std::vector<std::array<double,3>>& sphere, const double r_spher, const std::vector<Atom> atoms);
	
	const std::vector<std::array<double,3>>& get_points(){
		return points;
	}
	
	const size_t size(){
		return points.size();
	}
	
	friend std::ostream& operator<< (std::ostream& os, const grid& p);
	
	private:
	double vdwscl = 1.4;
	double vdwinc = 0.2;
	size_t nlayer = 6;
	
	std::vector<std::array<double,3>> points;
};
#endif
