#include "properties.h"

//using namespace libint2;

namespace libint2 {
int nthreads;

/// fires off \c nthreads instances of lambda in parallel
template <typename Lambda>
void parallel_do(Lambda& lambda) {
#ifdef _OPENMP
omp_set_num_threads(nthreads);
#pragma omp parallel
  {
    auto thread_id = omp_get_thread_num();
    lambda(thread_id);
  }
#else  // use C++11 threads
	std::vector<std::thread> threads;
	for (int thread_id = 0; thread_id != libint2::nthreads; ++thread_id) {
		threads.push_back(std::thread(lambda, thread_id));		
	}  // threads_id
	for (int thread_id = 0; thread_id < nthreads; ++thread_id)
		threads[thread_id].join();
		
	//std::cout << "work carried out on" << CPU_COUNT(&cpuset) << " cpus" << std::endl;
#endif
}
}


/*
void qd_calc (QM_residue& A){
		
	auto max_nprim = BasisSet::max_nprim(A.basis);
	auto max_l = BasisSet::max_l(A.basis);
		
	const std::vector<libint2::Shell>& obs = A.basis;	
	initialize();
	
	Engine s_engine(Operator::emultipole1, max_nprim, max_l);
	
	auto shell2bf = BasisSet::compute_shell2bf(A.basis);
	
	const auto& buf_vec = s_engine.results();
	const auto tstart = std::chrono::high_resolution_clock::now();
	
	double q{0}, dx{0}, dy{0}, dz{0};
	for(size_t s1=0; s1!=obs.size(); ++s1) {
		for(size_t s2=0; s2!=obs.size(); ++s2) {
		s_engine.compute(obs[s1], obs[s2]);
		auto s_shellset = buf_vec[0]; 
		auto x_shellset = buf_vec[1]; 
		auto y_shellset = buf_vec[2]; 
		auto z_shellset = buf_vec[3]; 
		if (s_shellset == nullptr){
			continue;  
		}
		auto bf1 = shell2bf[s1];  
		auto n1 = obs[s1].size(); 
		auto bf2 = shell2bf[s2];  
		auto n2 = obs[s2].size(); 
	
		// this iterates over integrals in this order
		for(size_t f1=0; f1!=n1; ++f1){
			for(size_t f2=0; f2!=n2; ++f2){				
				q += A.dm01(bf1 + f1, bf2 + f2)*s_shellset[f1*n2 + f2];
				dx += A.dm01(bf1 + f1, bf2 + f2)*x_shellset[f1*n2 + f2];
				dy += A.dm01(bf1 + f1, bf2 + f2)*y_shellset[f1*n2 + f2];
				dz += A.dm01(bf1 + f1, bf2 + f2)*z_shellset[f1*n2 + f2];
			}
		}
		}
	}
	
	const auto tstop = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> time_elapsed = tstop - tstart;

    std::cout << "Q = " << q << std::endl;
    std::cout << "d = {" << dx << ',' << dy << ',' << dz << '}' << std::endl;
    std::cout << "done (" << time_elapsed.count() << " s)" << std::endl;
	
	finalize();
}
*/

void v_calc (std::vector<libint2::Shell>& obs, Square_Matrix& dm){
	
	auto nthreads_cstr = getenv("LIBINT_NUM_THREADS");
    libint2::nthreads = 1;
    if (nthreads_cstr && strcmp(nthreads_cstr, "")) {
    std::istringstream iss(nthreads_cstr);
    iss >> libint2::nthreads;
    if (libint2::nthreads > 1 << 16 || libint2::nthreads <= 0) libint2::nthreads = 1;
	}
	
	size_t npoint{1000};
	
	std::random_device r;
	std::default_random_engine e1(10);
    std::uniform_real_distribution<double> uniform_distR(-100., 100.);
	
	std::vector<std::array<double, 3>> coords(npoint);
	for(auto& coord : coords){
		coord = {uniform_distR(e1),  uniform_distR(e1),  uniform_distR(e1)};
	}
		
	auto max_nprim = libint2::BasisSet::max_nprim(obs);
	auto max_l = libint2::BasisSet::max_l(obs);
				
	std::vector<libint2::Engine> engines(libint2::nthreads);
			
	for(auto& engine : engines){
		engine = libint2::Engine(libint2::Operator::nuclear, max_nprim, max_l);
	}
	
	auto shell2bf= libint2::BasisSet::compute_shell2bf(obs);
	
	auto v_lambda = [&](size_t thread_id) {
		double v;
		engines[thread_id] = libint2::Engine(libint2::Operator::nuclear, max_nprim, max_l);
		auto shell2bf = libint2::BasisSet::compute_shell2bf(obs);
		const auto& buf_vec = engines[thread_id].results();
			
		for(size_t w = 0; w < npoint; w++){
			if(!(w%100) && !thread_id) std::cout << w << " points done" << std::endl;
			if(w%libint2::nthreads != thread_id) continue;
											
			std::vector<std::pair<double, std::array<double, 3>>> p{
																		{1., 
																		{coords[w][0], coords[w][1], coords[w][2]}}
																	};
			engines[thread_id].set_params(p);
			
			v = 0.;
			for(size_t s1 = 0, s1max = obs.size(); s1 != s1max; s1++) {
				auto bf1 = shell2bf[s1];  
				auto n1 = obs[s1].size(); 
				for(size_t s2 = 0, s2max = obs.size(); s2 != s2max; s2++) {
					auto bf2 = shell2bf[s2];  
					auto n2 = obs[s2].size();
					 
					engines[thread_id].compute(obs[s1], obs[s2]);
					
					auto v_shellset = buf_vec[0]; 
					if (v_shellset == nullptr) continue;
					
					// this iterates over integrals in this order
					for(size_t f1 = 0; f1 != n1; f1++){
						for(size_t f2 = 0; f2 != n2; f2++){				
							v += dm(bf1 + f1, bf2 + f2)*v_shellset[f1*n2 + f2];
						}
					}
					
				}
			}
//			std::cout << v <<std::endl;
		}
	};
	
	const auto tstart = std::chrono::high_resolution_clock::now();
	
	libint2::parallel_do(v_lambda);
	
	const auto tstop = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> time_elapsed = tstop - tstart;
    
   
    
    std::cout << "done (" << time_elapsed.count() << " s)" << std::endl;
          
}
