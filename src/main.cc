#include <iomanip>

#include "qm_residue.h"
#include "properties.h"
#include "grid.h"
#include "fit.h"

#define P15 std::setprecision(15)

int main(int argc, char *argv[])
{
	std::string inpfile;
	std::string outfile;
	size_t nstate = 1;
	for(size_t i = 0; i < static_cast<size_t>(argc); ++i){
        if(std::string{argv[i]} == "-inp"){
			inpfile = std::string{argv[i+1]} + ".out";
			outfile = std::string{argv[i+1]} + ".esp";
		}
		if(std::string{argv[i]} == "-nst"){
			nstate = atol(argv[i+1]);
			if(nstate < 0) nstate = 1;
		}
    }
    
    std::cout << nstate << std::endl;
   
	// read QM file
	QM_residue p1(inpfile, nstate);	
	
// make grid
	polyhedron35 ico(3);
	grid grid0(ico.get_vertices(), ico.r_spher, p1.get_atoms());
/*		
	std::cout << "npoints = " << grid0.size() << std::endl;
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
	
// write output
	FILE* f = fopen(outfile.c_str(), "w");
	const auto& atoms = p1.get_atoms();
	
	fprintf(f, "everything is in atomic units\n\n\n"); 
	fprintf(f, "number of atoms = %zu\n\n", atoms.size());
	fprintf(f, "//label    nuclear charge    x    y    z\n#atomR\n"); 
	for(size_t i = 0, i_max = atoms.size(); i < i_max; ++i){
		 fprintf(f, "XX%11.6f%20.15f%20.15f%20.15f\n", static_cast<double>(atoms[i].atomic_number), atoms[i].r[0], atoms[i].r[1], atoms[i].r[2]);
	}
	fprintf(f, "#endatomR\n\n\n");
	fprintf(f, "number of charges = %zu\n\n//x    y    z    charge\n#charges\n", atoms.size());

	for(size_t i = 0, i_max = atoms.size(); i < i_max; ++i){		
		fprintf(f, "%20.15f%20.15f%20.15f%20.15f\n", atoms[i].r[0], atoms[i].r[1], atoms[i].r[2], esp[i]);
	}
	fprintf(f, "#endcharges\n"); 
	
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
