#include <iomanip>

#include "qm_residue.h"
#include "properties.h"
#include "grid.h"
#include "fit.h"
#include "algebra.h"

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

	std::cout << "npoints = " << grid0.size() << std::endl;
//libint 
	libint2::initialize();
	
	qm_prop engine(p1.get_basis(), p1.get_dms());
	
	auto qd_res = engine.qd_calc();
	
	std::vector<double> q(nstate);
	std::vector<std::array<double, 3>> d(nstate);
	auto l = engine.l_calc();
	auto p = engine.p_calc();
	
	for(size_t st = 0; st < nstate; ++st){
		q[st] = std::get<0>(qd_res[st]);
		d[st] = std::get<1>(qd_res[st]);
		std::cout << "transition 0->" << st+1 << ":\n";
		std::cout << "q = " << q[st] << std::endl;
		std::cout << "dipole:" << std::endl;
		std::cout << d[st][0] << ' ' << d[st][1] << ' ' << d[st][2] << ' ' << std::endl;
		std::cout << "transition angular momentum:" << std::endl;
		std::cout << l[st][0] << ' ' << l[st][1] << ' ' << l[st][2] << ' ' << std::endl;
		std::cout << "transition momentum:" << std::endl;
		std::cout << p[st][0] << ' ' << p[st][1] << ' ' << p[st][2] << ' ' << std::endl;
	}
	 
	auto v = engine.v_calc(grid0.get_points());

	libint2::finalize();

	fit esp_engine;
	std::vector<Eigen::VectorXd> esp(nstate);
	
	for(size_t st = 0; st < nstate; ++st){
		esp[st] = esp_engine(p1.get_atoms(), grid0.get_points(), v[st]);	
	}
	
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
	fprintf(f, "number of charges = %zu\n\n", atoms.size());
	fprintf(f, "//x    y    z    charge\n#charges\n");

	for(size_t st = 0; st < nstate; ++st){
		fprintf(f, "transition 0->%zu:\n", st+1);
		for(size_t i = 0, i_max = atoms.size(); i < i_max; ++i){		
			fprintf(f, "%20.15f%20.15f%20.15f%20.15f\n", atoms[i].r[0], atoms[i].r[1], atoms[i].r[2], esp[st][i]);
		}
	}
	fprintf(f, "#endcharges\n"); 
	
	fprintf(f, "\n\nsingle point properties:\n");
	for(size_t st = 0; st < nstate; ++st){
		q[st] = std::get<0>(qd_res[st]);
		d[st] = std::get<1>(qd_res[st]);
		fprintf(f, "transition 0->%zu:\n", st+1);
		fprintf(f, "q = %20.15f\n", q[st]);
		fprintf(f, "dipole:\n");
		fprintf(f, "%20.15f%20.15f%20.15f\n", d[st][0], d[st][1], d[st][2]);
		fprintf(f, "angular momentum:\n");
		fprintf(f, "%20.15f%20.15f%20.15f\n", l[st][0], l[st][1], l[st][2]);
		fprintf(f, "momentum:\n");
		fprintf(f, "%20.15f%20.15f%20.15f\n", p[st][0], p[st][1], p[st][2]);
		fprintf(f, "diff = %20.15f\n", 1 - dot(d[st], p[st])/norm(d[st])/norm(p[st]));
	}
	
	fprintf(f, "fitted properties\n"); 
	for(size_t st = 0; st < nstate; ++st){
		double q_fit = 0.;
		std::array<double, 3> d_fit {0., 0., 0.};
		for(size_t i = 0, i_max = atoms.size(); i < i_max; ++i){
			q_fit += esp[st][i];
			for(size_t k = 0; k < 3; ++k)
				d_fit[k] += atoms[i].r[k]*esp[st][i];
		}
		fprintf(f, "transition 0->%zu:\n", st+1);
		fprintf(f, "q = %20.15f\n", q_fit);
		fprintf(f, "%20.15f%20.15f%20.15f\n", d_fit[0], d_fit[1], d_fit[2]);
	}
	
	fclose(f);
	
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
