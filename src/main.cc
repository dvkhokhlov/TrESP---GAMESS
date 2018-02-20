#include "qm_residue.h"
#include "properties.h"
#include "grid.h"

int main()
{
//	QM_residue p1("tests/gms_chlb.out");	
	
	//libint2::initialize();
	
//	v_calc(p1.get_basis(), p1.get_dm());
	
	//libint2::finalize();
	
	polyhedron35 ico;
	
	ico.tesselate();
	ico.tesselate();
	ico.tesselate();
	ico.tesselate();
	ico.tesselate();
	ico.tesselate();
	
	//std::cout << ico;


	std::ofstream file;
	file.open("/home/daniilkh/ico.csv");
	
	for(auto& vertex_i : ico.vertices){
		for(size_t j = 0; j < 3; j++){
			file.width(10);
			file.precision(6);
			if(j < 2) file << vertex_i[j] << ",";
			else file << vertex_i[j];
		}
		file<<std::endl;
	}
}
