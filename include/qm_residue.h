#ifndef QM_RESIDUE_H_INCLUDED
#define QM_RESIDUE_H_INCLUDED

#include "libint2.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <memory>
#include <limits>
#include <ctgmath>

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

class QM_residue
{
	public:
	// constructors
	QM_residue() = delete;
	QM_residue(const std::string& _qm_fname);
	
	~QM_residue() = default;
	
	inline std::vector<libint2::Shell>& get_basis()
	{
		return basis;
	}
	
	inline Square_Matrix& get_dm()
	{
		return dm01;
	}	
	
	friend void qd_calc (QM_residue&);
//	friend double v_calc1 (QM_residue&);
//	friend double v_calc2 (QM_residue&);
	
	private:
	// add reading mode!
	bool pure = false;
	
	size_t natom;
	size_t ncgto;
	size_t nshell;
	
	std::string qm_fname;
		
	std::vector<libint2::Shell> basis;
	std::vector<libint2::Atom> atoms;
	
	Square_Matrix dm01;

// auxiliary functions/variables
	bool parsQ = false;
	std::string tmp;
	std::ifstream qm_file;
	void read_shell (const std::streampos&, const std::streampos&, size_t atom);	
	void read_pars ();
	void read_atoms();
	void read_basis ();
	void read_ecxprp ();
	void resort_dm ();
};

#endif
