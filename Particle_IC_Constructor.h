#ifndef Particle_IC_Constructor_H
#define Particle_IC_Constructor_H

#define DEBUG
#include<random>

#include<fstream>
#include"vector"
#include<string.h>
#include<sstream>

#include<iostream>


#include<time.h>
using namespace std;

/***gsl library***/
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>


extern double Interpolation( const int N, const double Table_x[], const double Table_y[], const double x );
extern double Mis_test( const int N, const double Table_x[], const double Table_y[], const double x );
extern int LoadTable( double *&Data, const char *FileName, const int NCol_Target, const int TCol[],
                   const bool AllocMem );

#define size_Models 1000
typedef struct Models_Input_Parameter{
   double* Models_Rho0;
   double* Models_R0;
   double* Models_MaxR;
   double ** Models_Center;
   double ** Models_BulkVel;
   double* Models_GasMFrac;
   int*    Models_MassProfNBin;

   int Models_num;
   vector<string> Models_Paras;
   vector<string> Models_Type;
   vector<string> Models_Profile;
   double* Models_Alpha;
   int*    Models_r_col;
   int*    Models_rho_col;
   bool* Models_truncation;
}MP;


class Particle_IC_Constructor
{
    public:
        Particle_IC_Constructor();
        virtual ~Particle_IC_Constructor();
        void init(string type,double al,double newton_g,double rho,double r,int nbin,double rmax,bool trunc_flag,double trunc_fac,int r_col,int rho_col,const char* Filename);        
        double set_vel(double x);
        double set_radius();
        double set_vel_test(double r);
        double vel_dispersion(double x);
        double specific_prob(double x,double vm);

        //Initialization through Table
        void initialize_mass_UNKNOWN(int row);
        void initialize_pot_UNKNOWN(int row);

        //Other type of models
        void initialize_mass_others();
        void initialize_pot_others();

        void initialize_prob_dens();
        
        //statistics
        double slope(double* a,double* b,int start,int fin);
        void smooth_all(double* x,int start,int fin);
        double set_mass( double x );
        double set_rho( double x );

        //Random Number
        double randomReal(double low,double high);
        
        
        
        MP params;
    protected:
        
    private:
        string model_type;
        double MassProf_Models( const double r );
        double integration_eng_base_Models(double eng);

        double prob_dens[size_Models];
        double int_prob_dens[size_Models];
        double psi[size_Models];
        double delta;
        double eng_min_Models;

        //statistics
        double ave(double* a,int start,int fin);
        double var_n(double* a,int start,int fin);
        double cor(double* x,double* y,int start,int fin);
        void mask(double* x,int start,int fin);
        void add_num(double* x,int start,int fin);

        //truncation
        bool Trunc_Flag;
        double Trunc_Fac;


};

#endif //Particle_IC_Constructor_H
