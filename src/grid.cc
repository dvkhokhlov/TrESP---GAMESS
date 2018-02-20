#include "grid.h"

void polyhedron35::gen_edges()
{
	edges.clear();
	edges.shrink_to_fit();
	for(size_t i = 0, size = vertices.size(); i < size; i++){
		for(size_t j = i + 1; j < size; j++){
			if(fabs(dist2(vertices[i], vertices[j]) - unit_len2) < 1E-6) 
				edges.push_back(edge{
										false, 0, std::array<size_t, 2>{i, j}
									});
		}
	}
}

// ugly brute-force
void polyhedron35::gen_faces()
{
	faces.clear();
	
	for(size_t i = 0, size = edges.size(); i < size; i++){
		for(size_t j = i + 1; j < size; j++){
			for(size_t k = j + 1; k < size; k++){
				if(edges[i].v[0] == edges[j].v[0] 
				&& (edges[i].v[1] == edges[k].v[0]
				&& edges[j].v[1] == edges[k].v[1]))
					faces.push_back(face{
											i, j, k
				//							std::array<size_t, 3>{edges[i].v[0], edges[k].v[0], edges[j].v[1]},
				//							std::array<size_t, 3>{i, j, k}
										});
			}
		}
	}
}

polyhedron35::polyhedron35()
{
	std::array<double,3> crd{0., 1., t};
	
	for(size_t i = 0; i < 3; i++)
		for(size_t k = 0; k < 2; k++)
			for(size_t l = 0; l < 2; l++){
				std::array<double,3> crd2{crd[0], powm1(crd[1], k), powm1(crd[2], l)};
				vertices.push_back(std::array<double,3> {crd2[i%3], crd2[(i+1)%3], crd2[(i+2)%3]});
	}

	gen_edges();
	gen_faces();
			
}

void polyhedron35::tesselate()
{
	auto vertex_num = faces.size();
	
	// dividing edges
	
	for(size_t i = 0, size = edges.size(); i < size; i++){
		auto& edge_current = edges[i];
		auto p = mid_v(vertices[edge_current.v[0]], vertices[edge_current.v[1]]);
		dscal(p, r/norm(p));
		
		vertices.push_back(p);
		
		edge_current.spawned_v = vertex_num;
		
		vertex_num++;
	}

	size_t edges_num = 0;
	std::vector<edge> edges_new;	
	// generating new ones
	for(size_t i = 0, size = faces.size(); i < size; i++){
		for(size_t j = 0; j < 3; j++){
			auto& edge_current = edges[faces[i].e[j]];
			if(edge_current.spawned_v){
				edges_new.push_back(edge{
										false, 0,  std::array<size_t, 2>{edge_current.v[0], edge_current.spawned_v}
										});
				edges_new.push_back(edge{
										false, 0, std::array<size_t, 2>{edge_current.v[1], edge_current.spawned_v}
										});
													
			}
		}
		
		for(size_t j = 0; j < 3; j++){
			auto& edge1 = edges[faces[i].e[j%3]];
			auto& edge2 = edges[faces[i].e[(j+1)%3]];
			edges_new.push_back(edge{
									false, 0,  std::array<size_t, 2>{edge1.spawned_v, edge2.spawned_v}
									});
		}
	}
	
	edges = edges_new;

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
	
	os << p.edges.size() << " edges" << std::endl;
	for(auto& edge_i : p.edges){
		os << edge_i.v[0] << '\t' << edge_i.v[1] << std::endl;
	}
	

	os << p.faces.size() << " faces:" << std::endl;
/*	os << "point1\tpoint2\tpoint3" << std::endl;
	for(auto& face_i : p.faces){
		for(auto& point_i : face_i.v)
			os << point_i << '\t';
		os <<  std::endl;
	}*/
	os << "edge1\tedge2\tedge3" << std::endl;
	for(auto& face_i : p.faces){
		for(auto& edge_i : face_i.e)
			os << '{' << p.edges[edge_i].v[0] << ','  << p.edges[edge_i].v[1] << '}'  << '\t';
		os <<  std::endl;
	}
	
	return os;
	
}

