#ifndef FIT_H_INCLUDED
#define FIT_H_INCLUDED

#include "algebra.h"
#include <vector>
#include <Eigen/Dense>

class fit{
	public:
	Eigen::VectorXd operator() (const std::vector<Atom>& atoms, const std::vector<std::array<double,3>>& points, const std::vector<double>& V){
	
	auto dim = atoms.size();
	
	A.resize(dim, dim);

	for(size_t i = 0; i < dim; i++){
        for(size_t j = i; j < dim; j++){
            A(i,j) = 0.;
            for(const auto& point : points){
                A(i,j) += 1.0/(dist(point, atoms[i].r)*dist(point, atoms[j].r));
            }
            A(j, i) = A (i, j);
        }
    }
	
	X.resize(dim);
	for(size_t i = 0; i < dim;i++){
		X[i]=0.;
        for(size_t h = 0, h_max = points.size(); h < h_max; h++){
            X[i] += V[h]/dist(points[h], atoms[i].r);
        }
    }
    
    return A.ldlt().solve(X);
    
	}
	
	private:
	Eigen::MatrixXd A;
	Eigen::VectorXd X;
	
};

#endif
