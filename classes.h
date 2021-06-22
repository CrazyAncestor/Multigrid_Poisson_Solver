#ifndef CLASSES_H
#define CLASSES_H

#include <stdlib.h>
#include <iostream>
#include <math.h>
#include<chrono>
#include "vector"

#include "macro.h"
#include "simulation_option.h"

#ifdef OMP_PARALLEL
#include<omp.h>
#endif

using namespace std;

class matrix
{ 
    int dim;
    int cells;
    double *value;
    double h;
    int did_x[3];
  
    public:
        matrix(int n,double hi);
        ~matrix();
        
        void display();
        void Error(const matrix &b);
  
        void   SOR_smoothing(const matrix& rho,int steps);
        void   SOR_Exact( const matrix &rho, int steps );
        double averaging(int idx);
        matrix Restriction();
        double insertion(int idx,int dim_in);
        matrix Interpolation(int idx);
        double Residual_error(const matrix& rho);
        matrix Residual(const matrix& rho);
        matrix Laplacian();
  
        void init_density();
        void init_potential();

        void reset();
  
        double get_h();
        double get_dim();
        double get_value( int idx );
        void   set_value( int idx, double val );
        void   add_value( int idx, double val );
        void   input_answer(int idx,double ans);
        void   set_boundary(double bc);

        matrix operator+(const matrix&);
        matrix operator-(const matrix&);

}; // CLASS : matrix


#endif // #ifndef CLASSES_H
