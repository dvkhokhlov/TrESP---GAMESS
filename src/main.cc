#include "qm_residue.h"
#include "properties.h"

int main()
{
	QM_residue p1("CBB601_libint.out");	
	
	libint2::initialize();
	
	v_calc(p1.get_basis(), p1.get_dm());
	
	libint2::finalize();
	
}
