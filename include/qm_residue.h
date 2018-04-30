#ifndef QM_RESIDUE_H_INCLUDED
#define QM_RESIDUE_H_INCLUDED

#include "libint2.hpp"

#include <Eigen/Dense>

#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <memory>
#include <limits>
#include <ctgmath>
#include "chemistry.h"
/*
struct Square_Matrix 
{
	Square_Matrix() = default;
	Square_Matrix(size_t _size) : size(_size), m(_size*_size) {}
	~Square_Matrix(){}
	
	void resize (size_t newsize){
		size = newsize;
		m.resize(size*size);
	}
	
	inline double& operator() (size_t i, size_t j)
	{
		return m[i*size + j];
	}
	
	size_t size;
	std::vector<double> m;
};
*/
class QM_residue
{
	public:
	// constructors
	QM_residue() = delete;
	QM_residue(const std::string& _qm_fname, size_t nstate=1);
	
	~QM_residue() = default;
	
// getters
	const std::vector<libint2::Shell>& get_basis();
	const std::vector<Atom>& get_atoms();
	Eigen::MatrixXd& get_dm(size_t);
	std::vector<Eigen::MatrixXd>& get_dms();
	
	private:
//
	size_t nstate;
	
// residue data	
	size_t natom;
	size_t ncgto;
	size_t nshell;
	std::vector<libint2::Shell> basis;
	std::vector<Atom> atoms;
//	Square_Matrix dm01;
	std::vector<Eigen::MatrixXd> dm;

	
// TODO: add reading mode!
	bool pure = false;
	
// strings/IO variables	
	bool parsQ = false;
	std::string tmp;
	std::string qm_fname;
	std::ifstream qm_file;
	
// auxiliary functions for IO
	void read_shell (const std::streampos&, const std::streampos&, size_t atom);	
	void read_pars ();
	void read_atoms();
	void read_basis ();
	void read_ecxprp (size_t);
	
// resort and renormalize density matrix
	void resort_dm (size_t);
};

#endif
