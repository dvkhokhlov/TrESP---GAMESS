#include "grid.h"

inline double powm1 (double a, int pow)
{
	return pow%2 ? -a : a;
}

/*
 * 
 * Constrution of unit shell
 * 
 */

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
	
	// placing points on middles of edges
	// numbers of new points connected with older ones are stored in neigh matrix
	std::vector<std::vector<size_t>> neigh(vertices.size());
	for(size_t i = 0, i_max = adj_mtx.size(); i < i_max; i++){
		for(size_t j = 0, j_max = adj_mtx[i].size(); j < j_max; j++){
			auto& edge2 = adj_mtx[i][j]; 
			
			auto p = mid_v(vertices[i], vertices[edge2]);
			dscal(p, r_spher/norm(p));
			vertices.push_back(p);
			
			auto point_num = vertices.size() - 1;
			neigh[i].push_back(point_num);			
			neigh[edge2].push_back(point_num);			
		}
	}

	// partial filling of new adjacency matrix
	// filled only neighbors of old points
	std::vector<std::vector<size_t>> adj_mtx_new(vertices.size());
	for(size_t i = 0, i_max = neigh.size(); i < i_max; i++){
		for(size_t j = 0, j_max = neigh[i].size(); j < j_max; j++){
			adj_mtx_new[i].push_back(neigh[i][j]); 
		}
	}
	
	// search of trangular faces
	// stored in triples matrix
	std::vector<std::array<size_t, 3>> triples;
	for(size_t i = 0, i_max = adj_mtx.size(); i < i_max; i++){
		for(size_t j = 0, j_max = adj_mtx[i].size(); j < j_max; j++){
			auto& edge2 = adj_mtx[i][j];
			
			for(size_t k = 0, k_max = adj_mtx[i].size(); k < k_max; k++){
				for(size_t l = 0, l_max = adj_mtx[edge2].size(); l < l_max; l++){
					if(adj_mtx[i][k] == adj_mtx[edge2][l]){
						triples.push_back({i, edge2, adj_mtx[i][k]});
					}
				}
			}
			
		}
	}
	
	// search of tesselated triangles included in older ones
	for(auto& triple : triples){
		std::vector<size_t> triples_intern;
		for(size_t i = 0; i < 2; i++){
			for(size_t j = 0, j_max = neigh[triple[i]].size(); j < j_max; j++){
				auto nmatch{0};
				auto candidate = neigh[triple[i]][j];
				
				for(size_t k = 0, k_max = neigh[triple[i+1]].size(); k < k_max; k++){
						if(candidate == neigh[triple[i+1]][k]){nmatch++;}
				}

				if(i == 0){
					for(size_t k = 0, k_max = neigh[triple[i+2]].size(); k < k_max; k++){
						if(candidate == neigh[triple[i+2]][k]){nmatch++;}
					}
				}
				
				if(nmatch == 1){
					triples_intern.push_back(candidate);
				}
			}
		}
	
		if(triples_intern.size() != 3){
			throw std::runtime_error ("polyhedron35.tesselate(): number of point in tesselated triangle != 3");
		}
			
		for(size_t i = 0; i < 3; i++){
			for(size_t j = 0; j < 3; j++){
				if(triples_intern[i] < triples_intern[j]) adj_mtx_new[triples_intern[i]].push_back(triples_intern[j]);
			}
		}
		
	}

	adj_mtx = adj_mtx_new;
}

const std::vector<std::array<double,3>>& polyhedron35::get_vertices()
{
	return vertices;
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
	for(auto& vertex_i : p.vertices){
		for(size_t j = 0; j < 3; j++){
			os.width(10);
			os.precision(6);
			if(j < 2) os << vertex_i[j] << ",";
			else os << vertex_i[j];
		}
		os<<std::endl;
	}

	return os;
}

/*
 * 
 * Constrution of full grid
 * 
 */

grid::grid(const std::vector<std::array<double,3>>& sphere, double r_spher)
{
	for(size_t i = 0; i < nlayer; i++){
		auto vdwr = 1.0;
		auto factor = vdwr*vdwscl/r_spher;
		for(auto& point : sphere){
			auto p = point;
			dscal(p, factor);
			points.push_back(p);
		}
		vdwscl += vdwinc;
	}
}

grid::grid(const std::vector<std::array<double,3>>& sphere, const double r_spher, const std::vector<Atom> atoms)
{
	// non-symmetrical sparce adjacency matrix
	// interatomic distance = Ri*(vdwscl + (nlayer - 1)*vdwinc) + Rj*(vdwscl)
	std::vector<std::vector<size_t>> adj_mtx(atoms.size());
	for(size_t i = 0, i_max = atoms.size(); i < i_max; i++){
		auto rmax_i = vdwr(atoms[i].atomic_number)*(vdwscl + (nlayer - 1)*vdwinc);
		
		for(size_t j = 0, j_max = atoms.size(); j < j_max; j++){
			auto maxdist = rmax_i  + vdwr(atoms[i].atomic_number)*vdwscl;
			if(dist2(atoms[i].r, atoms[j].r) < maxdist*maxdist){
				adj_mtx[i].push_back(j);
			}
		}
		
	}
	
	for(size_t i = 0, i_max = atoms.size(); i < i_max; i++){
		auto vdw_r = vdwr(atoms[i].atomic_number);
		auto vdwscl_local = vdwscl;
		
		for(size_t j = 0; j < nlayer; j++){
			auto factor = vdw_r*vdwscl_local/r_spher;
			
			for(auto& point : sphere){
				auto p = point;
				dscal(p, factor);
				p[0] += atoms[i].r[0];
				p[1] += atoms[i].r[1];
				p[2] += atoms[i].r[2];
				
				auto crossed = false;
				for(size_t k = 0, k_max = adj_mtx[i].size(); k < k_max; k++){
					auto& atom_k = atoms[adj_mtx[i][k]];
					if(dist2(p, atom_k.r) < vdwr(atom_k.atomic_number)*vdwr(atom_k.atomic_number)){
						crossed = true;
						break;
					}
				}
				
				if(!crossed) {points.push_back(p);}
			}
			
			vdwscl_local += vdwinc;
		}
	}
}

std::ostream& operator<< (std::ostream& os, const grid& g)
{
	for(auto& point : g.points){
		for(size_t j = 0; j < 3; j++){
			os.width(10);
			os.precision(6);
			if(j < 2) os << point[j] << ",";
			else os << point[j];
		}
		os<<std::endl;
	}

	return os;
}
