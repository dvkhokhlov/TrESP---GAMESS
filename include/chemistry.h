#ifndef CHEMISTRY_H_INCLUDED
#define CHEMISTRY_H_INCLUDED
#include <stdexcept>

struct Atom {
	int atomic_number;
	std::array<double,3> r;
};

// van der Waals radii in Bohrs // Batsanov2001
inline double vdwr (int atomic_number){
	double vdwr;
	switch(atomic_number){
		case 1: vdwr = 1.45; break; 
		case 5: vdwr = 1.80; break;
		case 6: vdwr = 1.70; break;
		case 7: vdwr = 1.60; break;
		case 8: vdwr = 1.55; break;
		case 9: vdwr = 1.50; break;
		case 12: vdwr = 2.20; break;
		case 14: vdwr = 2.10; break;
		case 15: vdwr = 1.95; break;
		case 16: vdwr = 1.80; break;
		case 17: vdwr = 1.80; break;
		case 20: vdwr = 2.40; break;
		case 26: vdwr = 2.05; break;
		case 30: vdwr = 2.10; break;
		default: 
			std::cout << "atomic_number = " << atomic_number << std::endl;
			throw std::runtime_error("vdwr(): unsupported atomic number");
	}
	
	return vdwr*1.88973;
}
#endif


