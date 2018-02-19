#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <memory>
#include <limits>
#include <ctgmath>
#include <array>

struct edge
{
	bool dividedQ;
	size_t v1, v2;
};

struct face
{
	size_t v1, v2, v3;
};

class polyhedron35 
{
	public:
	
	polyhedron35() : vertices()
	{
		std::vector<double> edges(30);
		std::vector<double> faces(20);
		
		std::array<double,3> crd{0., 1., t};
		
		for(size_t i = 0; i < 3; i++)
			for(size_t k = 0; k < 2; k++)
				for(size_t l = 0; l < 2; l++){
					std::array<double,3> crd2{crd[0], powm1(crd[1], k), powm1(crd[2], l)};
					vertices.push_back(std::array<double,3> {crd2[i%3], crd2[(i+1)%3], crd2[(i+2)%3]});
		}
				
	};
	
	void print()
	{
		for(auto& vertex_i : vertices)
			std::cout << vertex_i[0] << ' ' << vertex_i[1] << ' ' << vertex_i[2] << std::endl;
	}
	
	private:
	inline double powm1 (double a, int pow)
	{
		return pow%2 ? -a : a;
	}
	
	double t = (1. + sqrt(5.))/2.;
	std::vector<std::array<double,3>> vertices;
	std::vector<edge> edges;
	std::vector<face> faces;
};

#endif
