#ifndef PARTICLE_H
#define PARTICLE_H

#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "vector"
#include "macro.h"
#include "simulation_option.h"

#ifdef OMP_PARALLEL
#include <omp.h>
#endif

using namespace std;
class particle
{
   double par_mass;
   double *par_pos;
   double *par_vel;
   double *par_acc;

   public:
       particle();
       particle( double mass, double *pos, double *vel );
       ~particle();

       void Par_UpdateAll( const double *vel, const double *acc, const double dt );
       bool Par_InBox();

       void Par_AddMassToCell( matrix &source );
       void Par_AddMassToCell_NGP( matrix &source, const int *pos_idx, const double *dist_to_left );
       void Par_AddMassToCell_CIC( matrix &source, const int *pos_idx, const double *dist_to_left );
       void Par_AddMassToCell_TSC( matrix &source, const int *pos_idx, const double *dist_to_left );

       void Par_SumAcc( double *acc, double **force );
       void Par_SumAcc_NGP( double *acc, double **force, const int *pos_idx, const double *dist_to_left );
       void Par_SumAcc_CIC( double *acc, double **force, const int *pos_idx, const double *dist_to_left );
       void Par_SumAcc_TSC( double *acc, double **force, const int *pos_idx, const double *dist_to_left );

       void Par_SetMass( const double m );
       void Par_SetPos( const double *pos );
       void Par_SetVel( const double *vel );
       void Par_SetAcc( const double *acc );
       
       void Par_GetMass( double &m );
       void Par_GetPos( double *pos );
       void Par_GetVel( double *vel );
       void Par_GetAcc( double *acc );

}; // CLASS : particle


#endif // #ifndef PARTICLE_H
