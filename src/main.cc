#include "qm_residue.h"
#include "properties.h"
#include "grid.h"

int main()
{

	QM_residue p1("tests/gms_7amc.out");	
	
	polyhedron35 ico;
	ico.n_tesselate(3);
	
	grid grid0(ico.get_vertices(), ico.r_spher, p1.get_atoms());
	
	libint2::initialize();
	
	qd_calc(p1.get_basis(), p1.get_dm());
	
	auto v = v_calc(grid0.get_points(), p1.get_basis(), p1.get_dm());
	
	libint2::finalize();
	
	std::ofstream file;
	file.open("grid.txt");

	auto pts = grid0.get_points();
	for(size_t i = 0, i_max = grid0.size(); i < i_max; i++){
		file << pts[i][0] << ' ' << pts[i][1] << ' ' << pts[i][2] << ' ' << v[i] << std::endl;
	}
	
}
