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



std::vector<std::tuple<double, std::array<double,3>>> qm_prop::qd_calc ()
{
	auto max_nprim = libint2::BasisSet::max_nprim(obs);
	auto max_l = libint2::BasisSet::max_l(obs);
			
	libint2::Engine s_engine(libint2::Operator::emultipole1, max_nprim, max_l);
	
	auto shell2bf = libint2::BasisSet::compute_shell2bf(obs);
	
	const auto& buf_vec = s_engine.results();
	
	auto nstate = dms.size();
	
	std::vector<double> q(nstate), dx(nstate), dy(nstate), dz(nstate);
	
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
				
				for(size_t st = 0; st < nstate; ++st){				
					q[st] += dms[st](bf1 + f1, bf2 + f2)*s_shellset[f1*n2 + f2];
					dx[st] += dms[st](bf1 + f1, bf2 + f2)*x_shellset[f1*n2 + f2];
					dy[st] += dms[st](bf1 + f1, bf2 + f2)*y_shellset[f1*n2 + f2];
					dz[st] += dms[st](bf1 + f1, bf2 + f2)*z_shellset[f1*n2 + f2];
				}
			}
		}
		}
	}

	std::vector<std::tuple<double, std::array<double,3>>> result(nstate);
	for(size_t st = 0; st < nstate; ++st){
		result[st] = std::make_tuple(q[st], std::array<double, 3>{dx[st], dy[st], dz[st]});
	}
		
	return std::move(result);
}


std::vector<std::vector<double>> qm_prop::v_calc (const std::vector<std::array<double,3>>& points)
{
	auto max_nprim = libint2::BasisSet::max_nprim(obs);
	auto max_l = libint2::BasisSet::max_l(obs);
				
	std::vector<libint2::Engine> engines(nthreads);
			
	for(auto& engine : engines){
		engine = libint2::Engine(libint2::Operator::nuclear, max_nprim, max_l);
	}
	
	auto shell2bf= libint2::BasisSet::compute_shell2bf(obs);
	
	auto nstate = dms.size();
	
	std::vector<std::vector<double>> v(nstate);
	for(size_t st = 0; st < nstate; ++st){
		v[st].resize(points.size());
	}
	
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
			
			for(size_t st = 0; st < nstate; ++st){
				v[st][w] = 0.;
			}
				
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
							
							for(size_t st = 0; st < nstate; ++st){				
								v[st][w] += dms[st](bf1 + f1, bf2 + f2)*v_shellset[f1*n2 + f2];
							}
								
						}
					}
					
				}
			}
		}
	};
		
	parallel_do(v_lambda, nthreads);
    
    return std::move(v);
}

std::vector<std::array<double,3>> qm_prop::l_calc ()
{
	auto max_nprim = libint2::BasisSet::max_nprim(obs);
	auto max_l = libint2::BasisSet::max_l(obs);
			
	libint2::Engine engine(libint2::Operator::emultipole1, max_nprim, max_l + 1);
	
	auto shell2bf = libint2::BasisSet::compute_shell2bf(obs);
	
	const auto& buf_vec = engine.results();

	static std::vector<std::string> der_name 
		{"x", "y", "z"};

	static std::vector<std::string> p_name 
		{"x", "y", "z"};
		
	static std::vector<std::string> d_name 
		{"xx", "xy", "xz", "yy", "yz", "zz"};
		
	static std::vector<std::string> f_name 
		{"xxx", "xxy", "xxz", "xyy", "xyz", "xzz", "yyy", "yyz", "yzz", "zzz"};	
		
	static std::map<std::string, size_t> p_map 
		{{"x", 0}, {"y", 1}, {"z", 2}};
				
	static std::map<std::string, size_t> d_map 
		{{"xx", 0}, {"xy", 1}, {"xz", 2}, {"yy", 3}, {"yz", 4}, {"zz", 5}};
		
	static std::map<std::string, size_t> f_map 
		{{"xxx", 0}, {"xxy", 1}, {"xxz", 2}, {"xyy", 3}, {"xyz", 4}, {"xzz", 5}, 
		 {"yyy", 6}, {"yyz", 7}, {"yzz", 8}, {"zzz", 9}};
		 
		 
	auto phi_plus = [=](int am, int orb, char der)-> int {		
		std::string str;
		switch(am){
			case 1: 
				str = std::string{der};
				std::sort(str.begin(), str.end());
				return p_map[str];
			case 2:
				str = p_name[orb] + std::string{der};
				std::sort(str.begin(), str.end());
				return d_map[str];
			case 3:
				str = d_name[orb] + std::string{der};
				std::sort(str.begin(), str.end());
				return f_map[str];
			default:
				throw std::runtime_error ("l_calc(): unsupported am > 2\n");
		}
	};

	auto phi_minus = [=](int am, int orb, char der) ->  int{		
		std::string str;
		std::string orb_name;
		auto pos = std::string::npos;
		switch(am){
			case 1: 
				str = std::string{der};
				std::sort(str.begin(), str.end());
				return p_name[orb] == str ? 0 : -1;
			case 2:
				str = d_name[orb] + std::string{der};
				std::sort(str.begin(), str.end());
				orb_name = d_name[orb];
				pos = orb_name.find(der);
				if(pos != std::string::npos) orb_name.erase(pos, 1);
				return pos != std::string::npos ? p_map[orb_name] : -1;
			default:
				throw std::runtime_error ("l_calc(): unsupported am > 2\n");
		}
	};
	
	auto nstate = dms.size();
	
	std::vector<std::array<double,3>> l(nstate);
	
	for(size_t st = 0; st < nstate; ++st){
		l[st][0] = 0.;
		l[st][1] = 0.;
		l[st][2] = 0.;
	}
	
	for(size_t s1=0; s1!=obs.size(); ++s1) {
		for(size_t s2=0; s2!=obs.size(); ++s2) {
			
		auto aux_shell_plus = obs[s2];
		size_t am = obs[s2].contr[0].l;
		
		aux_shell_plus.contr[0].l += 1;
		size_t am_plus = aux_shell_plus.contr[0].l;
		
		for(size_t i = 0, i_max = aux_shell_plus.contr[0].coeff.size(); i < i_max; ++i)
			aux_shell_plus.contr[0].coeff[i] *= -2*aux_shell_plus.alpha[i];
			
		auto bf1 = shell2bf[s1];  
		auto n1 = obs[s1].size(); 
		auto bf2 = shell2bf[s2];  
		auto n2 = obs[s2].size(); 
			
		std::map<std::string, Eigen::MatrixXd> s1_s2 
			{{"xdy", Eigen::MatrixXd::Zero(n1, n2)}, 
		 	 {"xdz", Eigen::MatrixXd::Zero(n1, n2)}, 
		 	 {"ydx", Eigen::MatrixXd::Zero(n1, n2)}, 
		 	 {"ydz", Eigen::MatrixXd::Zero(n1, n2)}, 
		 	 {"zdx", Eigen::MatrixXd::Zero(n1, n2)}, 
		 	 {"zdy", Eigen::MatrixXd::Zero(n1, n2)}};

		engine.compute(obs[s1], aux_shell_plus);
		
		std::map<char, const double*> shset_map 
			{{'x', buf_vec[1]},
			 {'y', buf_vec[2]},
			 {'z', buf_vec[3]}};
			 
		for(size_t i = 0; i < n1; ++i){
			for(size_t j = 0; j < n2; ++j){
									
					for(auto& mtx : s1_s2){
						mtx.second(i, j) = shset_map[mtx.first[0]][i*aux_shell_plus.size() + phi_plus(am_plus, j, mtx.first[2])];
					}	
			}
		}
		
		if(obs[s2].contr[0].l > 0){
			auto aux_shell_minus = obs[s2];
			aux_shell_minus.contr[0].l -= 1;
			
			engine.compute(obs[s1], aux_shell_minus);
		
			shset_map = std::map<char, const double*>
				{{'x', buf_vec[1]},
				 {'y', buf_vec[2]},
				 {'z', buf_vec[3]}};
			
			for(size_t i = 0; i < n1; ++i){
				for(size_t j = 0; j < n2; ++j){
						
						for(auto& mtx : s1_s2){
							auto orb_pos = phi_minus(am, j, mtx.first[2]);
							if(orb_pos != -1){
								mtx.second(i, j) += shset_map[mtx.first[0]][i*aux_shell_minus.size() + orb_pos];
							}
						}
				}
			}
		}
			
		for(size_t f1=0; f1!=n1; ++f1){
			for(size_t f2=0; f2!=n2; ++f2){
				
				for(size_t st = 0; st < nstate; ++st){								
					l[st][0] += dms[st](bf1 + f1, bf2 + f2)*(s1_s2["ydz"](f1, f2) - s1_s2["zdy"](f1, f2));
					l[st][1] += dms[st](bf1 + f1, bf2 + f2)*(s1_s2["zdx"](f1, f2) - s1_s2["xdz"](f1, f2));
					l[st][2] += dms[st](bf1 + f1, bf2 + f2)*(s1_s2["xdy"](f1, f2) - s1_s2["ydx"](f1, f2));
				}
			}
		}
		
		
		}
	}
		
	return std::move(l);
}

std::vector<std::array<double,3>> qm_prop::p_calc ()
{
	auto max_nprim = libint2::BasisSet::max_nprim(obs);
	auto max_l = libint2::BasisSet::max_l(obs);
			
	libint2::Engine engine(libint2::Operator::overlap, max_nprim, max_l + 1);
	
	auto shell2bf = libint2::BasisSet::compute_shell2bf(obs);
	
	const auto& buf_vec = engine.results();

	static std::vector<std::string> der_name 
		{"x", "y", "z"};

	static std::vector<std::string> p_name 
		{"x", "y", "z"};
		
	static std::vector<std::string> d_name 
		{"xx", "xy", "xz", "yy", "yz", "zz"};
		
	static std::vector<std::string> f_name 
		{"xxx", "xxy", "xxz", "xyy", "xyz", "xzz", "yyy", "yyz", "yzz", "zzz"};	
		
	static std::map<std::string, size_t> p_map 
		{{"x", 0}, {"y", 1}, {"z", 2}};
				
	static std::map<std::string, size_t> d_map 
		{{"xx", 0}, {"xy", 1}, {"xz", 2}, {"yy", 3}, {"yz", 4}, {"zz", 5}};
		
	static std::map<std::string, size_t> f_map 
		{{"xxx", 0}, {"xxy", 1}, {"xxz", 2}, {"xyy", 3}, {"xyz", 4}, {"xzz", 5}, 
		 {"yyy", 6}, {"yyz", 7}, {"yzz", 8}, {"zzz", 9}};
		 
		 
	auto phi_plus = [=](int am, int orb, char der)-> int {		
		std::string str;
		switch(am){
			case 1: 
				str = std::string{der};
				std::sort(str.begin(), str.end());
				return p_map[str];
			case 2:
				str = p_name[orb] + std::string{der};
				std::sort(str.begin(), str.end());
				return d_map[str];
			case 3:
				str = d_name[orb] + std::string{der};
				std::sort(str.begin(), str.end());
				return f_map[str];
			default:
				throw std::runtime_error ("l_calc(): unsupported am > 2\n");
		}
	};

	auto phi_minus = [=](int am, int orb, char der) ->  int{		
		std::string str;
		std::string orb_name;
		auto pos = std::string::npos;
		switch(am){
			case 1: 
				str = std::string{der};
				std::sort(str.begin(), str.end());
				return p_name[orb] == str ? 0 : -1;
			case 2:
				str = d_name[orb] + std::string{der};
				std::sort(str.begin(), str.end());
				orb_name = d_name[orb];
				pos = orb_name.find(der);
				if(pos != std::string::npos) orb_name.erase(pos, 1);
				return pos != std::string::npos ? p_map[orb_name] : -1;
			default:
				throw std::runtime_error ("l_calc(): unsupported am > 2\n");
		}
	};

	auto nstate = dms.size();
	
	std::vector<std::array<double,3>> p(nstate);
	
	for(size_t st = 0; st < nstate; ++st){
		p[st][0] = 0.;
		p[st][1] = 0.;
		p[st][2] = 0.;
	}
	
	
	for(size_t s1=0; s1!=obs.size(); ++s1) {
		for(size_t s2=0; s2!=obs.size(); ++s2) {
			
		auto aux_shell_plus = obs[s2];
		size_t am = obs[s2].contr[0].l;
		
		aux_shell_plus.contr[0].l += 1;
		size_t am_plus = aux_shell_plus.contr[0].l;
		
		for(size_t i = 0, i_max = aux_shell_plus.contr[0].coeff.size(); i < i_max; ++i)
			aux_shell_plus.contr[0].coeff[i] *= -2*aux_shell_plus.alpha[i];
			
		auto bf1 = shell2bf[s1];  
		auto n1 = obs[s1].size(); 
		auto bf2 = shell2bf[s2];  
		auto n2 = obs[s2].size(); 
			
		std::map<std::string, Eigen::MatrixXd> s1_s2 
			{{"dx", Eigen::MatrixXd::Zero(n1, n2)},
			 {"dy", Eigen::MatrixXd::Zero(n1, n2)}, 
		 	 {"dz", Eigen::MatrixXd::Zero(n1, n2)}}; 

		engine.compute(obs[s1], aux_shell_plus);
		
		auto shset_map = buf_vec[0];
		
		for(size_t i = 0; i < n1; ++i){
			for(size_t j = 0; j < n2; ++j){
									
					for(auto& mtx : s1_s2){
						mtx.second(i, j) = shset_map[i*aux_shell_plus.size() + phi_plus(am_plus, j, mtx.first[1])];
					}	
			}
		}
		
		if(obs[s2].contr[0].l > 0){
			auto aux_shell_minus = obs[s2];
			aux_shell_minus.contr[0].l -= 1;
			
			engine.compute(obs[s1], aux_shell_minus);
		
			shset_map = buf_vec[0];
			
			for(size_t i = 0; i < n1; ++i){
				for(size_t j = 0; j < n2; ++j){
						
						for(auto& mtx : s1_s2){
							auto orb_pos = phi_minus(am, j, mtx.first[1]);
							if(orb_pos != -1){
								mtx.second(i, j) += shset_map[i*aux_shell_minus.size() + orb_pos];
							}
						}
				}
			}
		}
			
		for(size_t f1=0; f1!=n1; ++f1){
			for(size_t f2=0; f2!=n2; ++f2){
				for(size_t st = 0; st < nstate; ++st){								
					p[st][0] += dms[st](bf1 + f1, bf2 + f2)*s1_s2["dx"](f1, f2);
					p[st][1] += dms[st](bf1 + f1, bf2 + f2)*s1_s2["dy"](f1, f2);
					p[st][2] += dms[st](bf1 + f1, bf2 + f2)*s1_s2["dz"](f1, f2);
				}
			}
		}
		
		
		}
	}
		
	return std::move(p);
}
