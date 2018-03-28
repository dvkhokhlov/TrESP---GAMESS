#include "qm_residue.h"
#include "properties.h"
#include "grid.h"
#include "fit.h"

int main()
{
// read QM file
	QM_residue p1("tests/gms_7amc.out");	
	
// make grid
	polyhedron35 ico(3);
	grid grid0(ico.get_vertices(), ico.r_spher, p1.get_atoms());
		
//libint 
	libint2::initialize();
	
	qm_prop engine(p1.get_basis(), p1.get_dm());
	
	double q;
	std::array<double, 3> d;
	std::tie(q, d) = engine.qd_calc();
	auto v = engine.v_calc(grid0.get_points());
	
	libint2::finalize();
	
	fit esp_engine;
	
	auto esp = esp_engine(p1.get_atoms(), grid0.get_points(), v);
	
	for(size_t i = 0, i_max = esp.size(); i < i_max; i++)
		std::cout << esp[i] << std::endl;
		
// write
/*	
	std::ofstream file;
	file.open("grid.txt");

	auto pts = grid0.get_points();
	for(size_t i = 0, i_max = grid0.size(); i < i_max; i++){
		file << pts[i][0] << ' ' << pts[i][1] << ' ' << pts[i][2] << ' ' << v[i] << std::endl;
	}
	*/
}
