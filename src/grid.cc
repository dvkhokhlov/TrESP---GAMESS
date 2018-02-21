#include "grid.h"
inline double powm1 (double a, int pow)
{
	return pow%2 ? -a : a;
}

static double t = (1. + sqrt(5.))/2.;
static double unit_len2 = 4.;
static double r_spher = sqrt(1. + t*t);

polyhedron35::polyhedron35()
{
	std::array<double,3> crd{0., 1., t};
	
	for(size_t i = 0; i < 3; i++)
		for(size_t k = 0; k < 2; k++)
			for(size_t l = 0; l < 2; l++){
				std::array<double,3> crd2{crd[0], powm1(crd[1], k), powm1(crd[2], l)};
				vertices.push_back(std::array<double,3> {crd2[i%3], crd2[(i+1)%3], crd2[(i+2)%3]});
	}
	
	auto nvertex = vertices.size();
	adj_mtx.resize(nvertex);
	
	for(size_t i = 0, size = vertices.size(); i < size; i++){
		for(size_t j = i + 1; j < size; j++){
			if(fabs(dist2(vertices[i], vertices[j]) - unit_len2) < 1E-6)
				adj_mtx[i].push_back(j);
		}
	}	
}


void polyhedron35::tesselate()
{
	std::vector<std::vector<size_t>> neigh(vertices.size());
			
	for(size_t i = 0, i_max = adj_mtx.size(); i < i_max; i++){
		for(size_t j = 0, j_max = adj_mtx[i].size(); j < j_max; j++){
			auto& edge2 = adj_mtx[i][j]; 
			
			auto p = mid_v(vertices[i], vertices[edge2]);
			dscal(p, r_spher/norm(p));
			vertices.push_back(p);
						
			neigh[i].push_back(vertices.size() - 1);			
			neigh[edge2].push_back(vertices.size() - 1);			
		}
	}
/*	
	for(auto& row : neigh){
		for(auto& elem : row){
			std::cout << elem << '\t';
		}
		std::cout << std::endl;
	}
*/	
	std::vector<std::vector<size_t>> adj_mtx_new(vertices.size());
	
	for(size_t i = 0, i_max = neigh.size(); i < i_max; i++){
		for(size_t j = 0, j_max = neigh[i].size(); j < j_max; j++){
			adj_mtx_new[i].push_back(neigh[i][j]); 
		}
	}
	
	std::vector<std::array<size_t, 3>> triples;
	for(size_t i = 0, i_max = adj_mtx.size(); i < i_max; i++){
		for(size_t j = 0, j_max = adj_mtx[i].size(); j < j_max; j++){
			auto edge2 = adj_mtx[i][j];
			
			for(size_t k = 0, k_max = adj_mtx[i].size(); k < k_max; k++){
				for(size_t l = 0, l_max = adj_mtx[edge2].size(); l < l_max; l++){
					if(adj_mtx[i][k] == adj_mtx[edge2][l]){
//						std::cout << i << ' ' << edge2 << ' ' << adj_mtx[i][k] << std::endl;
						triples.push_back({i, edge2, adj_mtx[i][k]});
					}
				}	
			}
		}
	}
	
	std::vector<size_t> triples_intern;
	for(auto& triple : triples){
		std::vector<size_t> triples_intern;
		for(size_t i = 0; i < 2; i++){
			for(size_t k = 0, k_max = neigh[triple[i]].size(); k < k_max; k++){
				auto nmatch{0};
				auto candidate = neigh[triple[i]][k];
				
				for(size_t p = 0, p_max = neigh[triple[i+1]].size(); p < p_max; p++){
						if(candidate == neigh[triple[i+1]][p]){nmatch++;}
				}

				if(i == 0)
				for(size_t p = 0, p_max = neigh[triple[i+2]].size(); p < p_max; p++){
						if(candidate == neigh[triple[i+2]][p]){nmatch++;}
				}
				
				if(nmatch == 1){
					triples_intern.push_back(candidate);
				}
			}
		}
/*		
		if(triples_intern.size() != 3){
			std::cout << "error" << std::endl; getchar();
		}
*/		
		for(size_t i = 0; i < 3; i++){
			for(size_t j = 0; j < 3; j++){
				if(triples_intern[i] < triples_intern[j]) adj_mtx_new[triples_intern[i]].push_back(triples_intern[j]);
			}
		}
		
	}


	adj_mtx = adj_mtx_new;
}

void polyhedron35::n_tesselate(size_t n)
{
	assert(n < 11ul);
	
	for(size_t i = 0; i < n; i++)
		tesselate();
}

void polyhedron35::tesselate_to(size_t n)
{
	assert(n < 10000000ul);
	while(vertices.size() < n)
		tesselate();
}


std::ostream& operator<< (std::ostream& os, const polyhedron35& p)
{
	os << p.vertices.size() << " vertices" << std::endl;
	for(auto& vertex_i : p.vertices){
		for(size_t j = 0; j < 3; j++){
			os.width(10);
			os.precision(6);
			if(j < 2) os << vertex_i[j] << ",";
			else os << vertex_i[j];
		}
		os<<std::endl;
	}
/*	
	os << " adjacency matrix" << std::endl;
	for(size_t i = 0; i < p.adj_mtx.size(); i++){
		os << i << '\t';
		for(auto& elem : p.adj_mtx[i]){
			os << elem << '\t';
		}
		os << std::endl;
	}
*/	
	return os;
	
}

