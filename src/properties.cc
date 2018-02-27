#include "properties.h"

/// fires off \c nthreads instances of lambda in parallel
template <typename Lambda>
void parallel_do(Lambda& lambda, int nthreads) {
	std::vector<std::thread> threads;
	
	for (int thread_id = 0; thread_id != nthreads; ++thread_id) {
		threads.push_back(std::thread(lambda, thread_id));		
	}  // threads_id
	
	for (int thread_id = 0; thread_id < nthreads; ++thread_id)
		threads[thread_id].join();
}



std::tuple<double, std::array<double,3>> qm_prop::qd_calc ()
{
	auto max_nprim = libint2::BasisSet::max_nprim(obs);
	auto max_l = libint2::BasisSet::max_l(obs);
			
	libint2::Engine s_engine(libint2::Operator::emultipole1, max_nprim, max_l);
	
	auto shell2bf = libint2::BasisSet::compute_shell2bf(obs);
	
	const auto& buf_vec = s_engine.results();
	
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
				q += dm(bf1 + f1, bf2 + f2)*s_shellset[f1*n2 + f2];
				dx += dm(bf1 + f1, bf2 + f2)*x_shellset[f1*n2 + f2];
				dy += dm(bf1 + f1, bf2 + f2)*y_shellset[f1*n2 + f2];
				dz += dm(bf1 + f1, bf2 + f2)*z_shellset[f1*n2 + f2];
			}
		}
		}
	}
		
	return std::make_tuple(q, std::array<double, 3>{dx, dy, dz});
}


std::vector<double> qm_prop::v_calc (const std::vector<std::array<double,3>> points)
{
	auto max_nprim = libint2::BasisSet::max_nprim(obs);
	auto max_l = libint2::BasisSet::max_l(obs);
				
	std::vector<libint2::Engine> engines(nthreads);
			
	for(auto& engine : engines){
		engine = libint2::Engine(libint2::Operator::nuclear, max_nprim, max_l);
	}
	
	auto shell2bf= libint2::BasisSet::compute_shell2bf(obs);
	
	std::vector<double> v(points.size());
	
	auto v_lambda = [&](size_t thread_id) {
		engines[thread_id] = libint2::Engine(libint2::Operator::nuclear, max_nprim, max_l);
		auto shell2bf = libint2::BasisSet::compute_shell2bf(obs);
		const auto& buf_vec = engines[thread_id].results();
			
		for(size_t w = 0, w_max = points.size(); w < w_max; w++){
			if(!(w%100) && !thread_id) std::cout << w << " points done" << std::endl;
			if(w%nthreads != thread_id) continue;
											
			std::vector<std::pair<double, std::array<double, 3>>> p{
																		{1., 
																		{points[w][0], points[w][1], points[w][2]}}
																	};
			engines[thread_id].set_params(p);
			
			v[w] = 0.;
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
							v[w] += dm(bf1 + f1, bf2 + f2)*v_shellset[f1*n2 + f2];
						}
					}
					
				}
			}
		}
	};
		
	parallel_do(v_lambda, nthreads);
    
    return std::move(v);
}
