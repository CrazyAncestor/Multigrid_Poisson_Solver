#ifndef INITIALIZE_H
#define INITIALIZE_H


#include "macro.h"
#include "simulation_option.h"


bool Init_matrix( matrix &mat, particle *pars );

bool Init_SinWave( matrix &mat, particle *pars );

bool Init_TwoBody( matrix &mat, particle * pars );

bool Init_NBody( matrix &mat, particle * pars );

bool Init_SolarOrbit(matrix& mat, particle* pars);

#endif // #ifndef INITIALIZE_H
