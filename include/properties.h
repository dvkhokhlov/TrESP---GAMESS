#ifndef PROPERTIES_H_INCLUDED
#define PROPERTIES_H_INCLUDED
#include "omp.h"
#include "qm_residue.h"

#include <pthread.h>
#include <random>
#include <chrono>
#include <algorithm>
#include <cstdio>
#include <tuple>
#include <thread>

class qm_prop{
	public:
//	qm_prop(const std::vector<libint2::Shell>& _obs, Square_Matrix& _dm) : obs(_obs), dm(_dm)
	qm_prop(const std::vector<libint2::Shell>& _obs, std::vector<Eigen::MatrixXd>& _dms) : obs(_obs), dms(_dms)
	{
		auto nthreads_cstr = getenv("LIBINT_NUM_THREADS");
		nthreads = 1;
		if (nthreads_cstr && strcmp(nthreads_cstr, "")) {
		std::istringstream iss(nthreads_cstr);
		iss >> nthreads;
		if (nthreads > 1 << 16 || nthreads <= 0) nthreads = 1;
		}
	}
	
	std::vector<std::tuple<double, std::array<double,3>>> qd_calc ();
	std::vector<std::vector<double>> v_calc (const std::vector<std::array<double,3>>& points);
	std::vector<std::array<double,3>> l_calc ();
	std::vector<std::array<double,3>> p_calc ();
	
	private:
	int nthreads;
	const std::vector<libint2::Shell>& obs;
//	Square_Matrix& dm;
	std::vector<Eigen::MatrixXd>& dms;
};

#endif
