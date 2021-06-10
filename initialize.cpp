#include "classes.h"
#include "particle.h"

#include "initialize.h"



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
        const int i = idx%BOX_N;
        const int j = ( idx%(BOX_N*BOX_N) ) / BOX_N;
        const int k = idx/(BOX_N*BOX_N);
        
        const double x = BOX_DX*i;
        const double y = BOX_DX*j;
        const double z = BOX_DX*k;

        mat.set_value( idx, -2.*sin(x)*sin(y)*sin(z) );
     
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
bool Init_TwoBody( matrix &mat, particle *pars )
{
    if ( N_PARS != 2 )
    {
        printf("Number of particles need to be 2 in this simulation problem.\n");
        return false;
    }

    mat.reset();

    // particles
    double *pos = new double[N_DIMS];
    double *vel = new double[N_DIMS];
    for ( int p = 0; p < N_PARS; p++ )
    {
        if (p == 0)
        {
            pos[0] = BOX_L/2.0 + 1.0;
            pos[1] = BOX_L/2.0;
            vel[0] = 0.0;
            vel[1] = 0.5;
        } else
        {
            pos[0] = BOX_L/2.0 - 1.0;
            pos[1] = BOX_L/2.0;
            vel[0] = 0.0;
            vel[1] = -0.5;
        }
        #if ( N_DIMS == 3 )
        pos[2] = 0.0;
        vel[2] = 0.0;
        #endif

        pars[p].Par_SetMass( 1.0 );
        pars[p].Par_SetPos( pos );
        pars[p].Par_SetVel( vel );
        pars[p].Par_AddMassToCell( mat );

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

    // particles
    double *pos = new double[N_DIMS];
    double *vel = new double[N_DIMS];
    double d_temp = BOX_L/100.;
    for ( int p = 0; p < N_PARS; p++ )
    {
        pos[0] = ( p/100 ) * d_temp;
        pos[1] = ( p%100 ) * d_temp;
        vel[0] = 0.0;
        vel[1] = 0.0;

        #if ( N_DIMS == 3 )
        pos[2] = BOX_L/2.0;
        vel[2] = 0.0;
        #endif 

        pars[p].Par_SetMass( 1.0 );
        pars[p].Par_SetPos( pos );
        pars[p].Par_SetVel( vel );
        pars[p].Par_AddMassToCell( mat );

    } // for ( int p = 0; p < N_PARS; p++ )

    delete[] pos;
    delete[] vel;

    return true;
} /// FUNCTION : Init_NBody
