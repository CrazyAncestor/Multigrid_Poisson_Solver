#ifndef SOLVER_H
#define SOLVER_H

#include "macro.h"
#include "simulation_option.h"

matrix Solver_Potential( matrix phi, matrix dens );

void Solver_Force( matrix pot, double **force );

matrix SOR_Method( matrix phi, matrix dens );

matrix V_Cycle( matrix phi, matrix dens );

matrix W_Cycle( matrix phi, matrix dens);

matrix FAS_Method( matrix phi, matrix dens );

matrix FMG_Method( matrix phi, matrix dens, int n_fmg);

#endif // #ifndef SOLVER_H
