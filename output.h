#ifndef OUTPUT_H
#define OUTPUT_H

#include "macro.h"
#include "simulation_option.h"

void Output_matrix( matrix mat, const char filename[] );

void Output_particles( particle *pars, const char filename[] );

void Output_parameter();

#endif
