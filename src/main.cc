#include "qm_residue.h"
#include "properties.h"
#include "grid.h"

int main()
{
	QM_residue p1("tests/gms_7amc.out");	
	
	libint2::initialize();
	
	qd_calc(p1.get_basis(), p1.get_dm());
	v_calc(p1.get_basis(), p1.get_dm());
	
	libint2::finalize();

/*	
	polyhedron35 ico;
	
	
	ico.n_tesselate(4);
	
	//std::cout << ico;


	std::ofstream file;
	file.open("ico.csv");

	file << ico;
	*/
}
