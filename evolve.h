#ifndef EVOLVE_H
#define EVOLVE_H

#include "macro.h"
#include "simulation_option.h"

double Evolve_GetDt( particle *pars, const double t_now, const int counter, bool &out_stat, bool &end_stat );

void Evolve_UpdateParticle( matrix &dens, particle *pars, const double dt );

void Evolve_Euler( matrix &dens, particle *pars, const double dt );

void Evolve_KDK( matrix &dens, particle *pars, const double dt );

void Evolve_DKD( matrix &dens, particle *pars, const double dt );

void Evolve_RK4( matrix &dens, particle *pars, const double dt );

void Evolve_DKDOrbit(matrix& dens, particle* pars, const double dt);

#endif // #ifndef EVOLVE_H
