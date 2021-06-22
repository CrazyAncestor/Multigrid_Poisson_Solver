#include "classes.h"
#include "particle.h"

#include "initialize.h"

#include "Particle_IC_Constructor.h"

//--------------------------------------------------------------------------------
// Function    : Init_matrix
// Description : Initialize the matrix by different PROB_NUM
// Note        :
// Input       : mat  : matrix
//               pars : particles
// Output      : if the initializaion is done or not.
//--------------------------------------------------------------------------------
bool Init_matrix( matrix &mat, particle *pars )
{
    bool init_stat = false;
    #if ( PROB_NUM == PROB_SINWAVE )
    init_stat = Init_SinWave( mat, pars );
    #elif ( PROB_NUM == PROB_TWOBODY )
    init_stat = Init_TwoBody( mat, pars );
    #elif ( PROB_NUM == PROB_NBODY )
    init_stat = Init_NBody( mat, pars );
    #elif ( PROB_NUM == PROB_SOLARORBIT )
    init_stat = Init_SolarOrbit(mat, pars);
    #else
    printf("No such a problem number %d.\n", PROB_NUM);
    #endif

    return init_stat;

} // FUNCTION : Init_matrix( matrix &mat )



//--------------------------------------------------------------------------------
// Function    : Init_SinWave
// Description : Matrix value set as f(x, y) = -2 * sin(x) * sin(y) 
// Note        :
// Input       : mat : matrix
// Output      : if the initializaion is done or not.
//--------------------------------------------------------------------------------
bool Init_SinWave( matrix &mat, particle *pars )
{
    for( int idx = 0; idx < N_CELLS; idx++ ) 
    {
        
    #if ( N_DIMS == 2 )
        const int i = idx / BOX_N;
        const int j = idx % BOX_N;

        const double x = BOX_DX * i;
        const double y = BOX_DX * j;
        mat.set_value(idx, -2. * sin(x) * sin(y));
    #endif
    #if ( N_DIMS == 3 )
        const int dim_2 = BOX_N * BOX_N;
        const int i = idx / dim_2;
        const int j = (idx % dim_2) / BOX_N;
        const int k = (idx % dim_2) % BOX_N;

        const double x = BOX_DX*i;
        const double y = BOX_DX*j;
        const double z = BOX_DX*k;

        mat.set_value( idx, -3.*sin(x)*sin(y)*sin(z) );
    #endif
     } // for( int idx = 0; idx < N_CELLS; idx++ ) 
    
    return true;
} // FUNCTION : Init_SinWave



//--------------------------------------------------------------------------------
// Function    : Init_TwoBody
// Description : Problem number 2
// Note        :
// Input       : mat : matrix
// Output      : if the initializaion is done or not.
//--------------------------------------------------------------------------------
bool Init_TwoBody(matrix& mat, particle* pars)
{
    if (N_PARS != 2)
    {
        printf("Number of particles need to be 2 in this simulation problem.\n");
        return false;
    }

    mat.reset();

    // particles
    double* pos = new double[N_DIMS];
    double* vel = new double[N_DIMS];
    for (int p = 0; p < N_PARS; p++)
    {
        if (p == 0)
        {
            pos[0] = BOX_L / 2.0;
            pos[1] = BOX_L / 2.0;
            vel[0] = 0.0;
            vel[1] = 0.0;
        }
        else
        {
            pos[0] = BOX_L / 2.0 + 1.0;
            pos[1] = BOX_L / 2.0;
            vel[0] = 0.0;
            vel[1] = pow(2 * PI, -0.5);
        }
#if ( N_DIMS == 3 )
        pos[2] = BOX_L / 2.0;
        vel[2] = 0.0;
#endif

        pars[p].Par_SetMass(1.0);
        if (p == 1)pars[p].Par_SetMass(0.0);
        pars[p].Par_SetPos(pos);
        pars[p].Par_SetVel(vel);
        pars[p].Par_AddMassToCell(mat);

    } // for ( int p = 0; p < N_PARS; p++ )

    delete[] pos;
    delete[] vel;

    return true;
} /// FUNCTION : Init_TwoBody

//--------------------------------------------------------------------------------
// Function    : Init_NBody
// Description : Problem number 3
// Note        :
// Input       : mat : matrix
// Output      : if the initializaion is done or not.
//--------------------------------------------------------------------------------
bool Init_NBody( matrix &mat, particle *pars )
{
    if ( N_PARS != 10000 )
    {
        printf("Number of particles need to be 10000 in this simulation problem.\n");
        return false;
    }
    
    mat.reset();
    // Initialize Particle IC Constructor
    double Newton_G      =  1.0;    //Gravitational Constant
    double Rho0          =  1.0;     // peak density [unit: M mass of sun/kpc^3]
    double R0            =  1.0 ;                 // scale radius [unit: kpc]
    double MaxR          =  2.0;                    // maximum radius for particles [unit: kpc]
    int MassProfNBin   = 1000    ;             // number of radial bins in the mass profile table [1000]
    double Alpha         =1     ;               // alpha parameter for Eiasto model [1]
    int r_col         = 0         ;           // number of the column of the radius of density profile, when model is "UNKNOWN"  [0]
    int rho_col        =1          ;          // number of the column of the density of density profile, when model is "UNKNOWN"  [1]
    double truncation    = 0       ;             // whether to turn on a smoothy truncation function of density near MaxR [0]

    Particle_IC_Constructor constructor_Models;
    constructor_Models.init("Plummer",Alpha,Newton_G,Rho0,R0,MassProfNBin,MaxR,truncation,0.7,r_col,rho_col,"NONE");


    // particles
    double *pos = new double[N_DIMS];
    double *vel = new double[N_DIMS];
    double d_temp = BOX_L/100.;

    //compute mass
    double par_m =constructor_Models.set_mass( MaxR/R0 )/N_PARS;
    for ( int p = 0; p < N_PARS; p++ )
    {
        #if ( N_DIMS == 2 )
        pos[0] = ( p/100 ) * d_temp;
        pos[1] = ( p%100 ) * d_temp;
        vel[0] = 0.0;
        vel[1] = 0.0;
        #endif 

        #if ( N_DIMS == 3 )
        double r = constructor_Models.set_radius();
        double v = constructor_Models.set_vel(r)/pow(4*PI,0.5);
        
        //random direction
        double theta1 = constructor_Models.randomReal(0,PI);
        double phi1   = constructor_Models.randomReal(0,2*PI);

        double theta2 = constructor_Models.randomReal(0,PI);
        double phi2   = constructor_Models.randomReal(0,2*PI);

        pos[0] = r*sin(theta1)*cos(phi1)+BOX_L/2.0;
        pos[1] = r*sin(theta1)*sin(phi1)+BOX_L/2.0;
        pos[2] = r*cos(theta1)+BOX_L/2.0;

        vel[0] = v*sin(theta2)*cos(phi2);
        vel[1] = v*sin(theta2)*sin(phi2);
        vel[2] = v*cos(theta2);
        #endif 

        pars[p].Par_SetMass( par_m );
        pars[p].Par_SetPos( pos );
        pars[p].Par_SetVel( vel );
        pars[p].Par_AddMassToCell( mat );
    } // for ( int p = 0; p < N_PARS; p++ )

    delete[] pos;
    delete[] vel;

    return true;
} /// FUNCTION : Init_NBody



//--------------------------------------------------------------------------------
// Function    : Init_SolarOrbit
// Description : Problem number 4
// Note        :
// Input       : mat : matrix
// Output      : if the initializaion is done or not.
//--------------------------------------------------------------------------------
bool Init_SolarOrbit(matrix& mat, particle* pars)
{
    if (N_PARS != 2)
    {
        printf("Number of particles need to be 2 in this simulation problem.\n");
        return false;
    }

    mat.reset();

    // particles
    double* pos = new double[N_DIMS];
    double* vel = new double[N_DIMS];
    for (int p = 0; p < N_PARS; p++)
    {
        if (p == 0)
        {
            pos[0] = BOX_L / 2.0;
            pos[1] = BOX_L / 2.0;
            vel[0] = 0.0;
            vel[1] = 0;
            pars[p].Par_SetMass(1E+3);
        }
        else
        {
            pos[0] = BOX_L / 2.0 + 1;
            pos[1] = BOX_L / 2.0;
            vel[0] = 0;
            vel[1] = pow(2 * PI, 1);//for 2D  /*pow(2 * pi, -0.5) for 3D;*/
            pars[p].Par_SetMass(1);
        }
#if ( N_DIMS == 3 )
        pos[2] = BOX_L / 2.0;
        vel[2] = 0.0;
#endif


        pars[p].Par_SetPos(pos);
        pars[p].Par_SetVel(vel);
        pars[p].Par_AddMassToCell(mat);



    } // for ( int p = 0; p < N_PARS; p++ )

    delete[] pos;
    delete[] vel;

    return true;
} /// FUNCTION : Init_SolarOrbit