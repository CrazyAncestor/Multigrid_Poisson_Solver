#include "classes.h"
#include "particle.h"

#include "solver.h"
#include "evolve.h"


//--------------------------------------------------------------------------------
// Function    : Evolve_GetDt
// Description : 
// Note        :
// Input       : 
// Output      : 
//--------------------------------------------------------------------------------
double Evolve_GetDt( particle *pars, const double t_now, const int counter, bool &out_stat, bool &end_stat )
{
    double min_dt = 0.01;
    double out_time = counter * OUT_DT;

    if ( t_now + min_dt >= END_TIME )
    {
        min_dt = END_TIME - t_now;
        out_stat = true;
        end_stat = true;
    } else if ( t_now + min_dt >= out_time )
    {
        min_dt = out_time - t_now;
        out_stat = true;
        end_stat = false;
    } else 
    {
        out_stat = false;
        end_stat = false;
    }
    
    return min_dt;
}



//--------------------------------------------------------------------------------
// Function    : Evolve_UpdateParticle
// Description : 
// Note        :
// Input       : 
// Output      : 
//--------------------------------------------------------------------------------
void Evolve_UpdateParticle( matrix &dens, particle *pars, const double dt )
{
    #if ( EVOLVE == EULER )
    Evolve_Euler( dens, pars, dt ); 
    #elif ( EVOLVE == KDK )
    Evolve_KDK( dens, pars, dt ); 
    #elif ( EVOLVE == DKD )
    Evolve_DKD( dens, pars, dt ); 
    #elif ( EVOLVE == RK4 )
    Evolve_RK4( dens, pars, dt ); 
    #endif



} // FUNCTION : Evolve_UpdateParticle



//--------------------------------------------------------------------------------
// Function    : Evolve_Euler
// Description : 
// Note        :
// Input       : 
// Output      : 
//--------------------------------------------------------------------------------
void Evolve_Euler( matrix &dens, particle *pars, const double dt ) 
{
    double **force;
    force = new double *[N_DIMS];
    for ( int d = 0; d < N_DIMS; d++ )    force[d] = new double [N_CELLS];

    // 1. solve acc
    const double dx = BOX_DX;
    matrix pot( BOX_N, dx );
    pot.init_potential();
    matrix solved_pot = Solver_Potential( pot, dens );
    Solver_Force( solved_pot, force );

    // 2. update
    double *vel = new double[N_DIMS];
    double *acc = new double[N_DIMS];
    
    for ( int p = 0; p < N_PARS; p++ )
    {
        pars[p].Par_GetVel( vel );
        pars[p].Par_SumAcc( acc, force ); // get the velocity
        //printf("P = %d, ax = %.5e ,ay = %.5e", p, acc[0], acc[1]);
        pars[p].Par_UpdateAll( vel, acc, dt );
    } // for ( int p = 0; p < N_PARS; p++ )
    //printf("\n");

    delete[] vel;
    delete[] acc;
    
    for ( int d = 0; d < N_DIMS; d++ )    delete[] force[d];
    delete[] force;

} // FUNCTION : Evolve_Euler



//--------------------------------------------------------------------------------
// Function    : Evolve_KDK
// Description : 
// Note        :
// Input       : 
// Output      : 
//--------------------------------------------------------------------------------
void Evolve_KDK( matrix &dens, particle *pars, const double dt )
{
} // FUNCTION : Evolve_KDK





//--------------------------------------------------------------------------------
// Function    : Evolve_DKD
// Description : 
// Note        : uncommend first when start to calculate orbit(156-158,190-192)
// Input       : 
// Output      : 
//--------------------------------------------------------------------------------
void Evolve_DKD(matrix& dens, particle* pars, const double dt)
{
    double** force;
    force = new double* [N_DIMS];
    for (int d = 0; d < N_DIMS; d++)    force[d] = new double[N_CELLS];

    double* pos = new double[N_DIMS];
    double* vel = new double[N_DIMS];
    double* acc = new double[N_DIMS];
    double r;
    double relativity_factor;
    int alpha;
    const double dx = BOX_DX;

    alpha = 0.0008;

    // 1. kick
    for (int p = 0; p < N_PARS; p++)
    {
        pars[p].Par_GetPos(pos);
        pars[p].Par_GetVel(vel);
        //printf("P = %d, px = %.5e ,py = %.5e ", p, pos[0], pos[1]);
        for (int d = 0; d < N_DIMS; d++)
        {
            //open it only when calculating orbit
            /*r = pow((pos[d] * pos[d] + pos[d + 1] * pos[d + 1]), 0.5);
            relativity_factor = 1 + alpha / (pow(r, 2));
            vel[d] -= (4 * PI * PI * pos[d] * 0.5 * dt) * relativity_factor / (pow(r, 3));*/
            pos[d] += vel[d] * 0.5 * dt;
        }
        //printf("=> px = %.5e ,py = %.5e ", pos[0], pos[1]);
        pars[p].Par_SetPos(pos);
    } // for ( int p = 0; p < N_PARS; p++ )
    //printf("\n");


    // 2. drift
    matrix pot(BOX_N, dx);
    pot.init_potential();
    matrix solved_pot = Solver_Potential(pot, dens);
    Solver_Force(solved_pot, force);

    for (int p = 0; p < N_PARS; p++)
    {
        pars[p].Par_SumAcc(acc, force); // get the velocity
        pars[p].Par_GetVel(vel);
        for (int d = 0; d < N_DIMS; d++)    vel[d] += acc[d] * dt;
        pars[p].Par_SetVel(vel);
    } // for ( int p = 0; p < N_PARS; p++ )


    // 3. kick
    for (int p = 0; p < N_PARS; p++)
    {
        pars[p].Par_GetPos(pos);
        pars[p].Par_GetVel(vel);
        for (int d = 0; d < N_DIMS; d++)    
        {
            //open it only when calculating orbit
            /*r = pow((pos[d] * pos[d] + pos[d + 1] * pos[d + 1]), 0.5);
            relativity_factor = 1 + alpha / (pow(r, 2));
            vel[d] -= (4 * PI * PI * pos[d] * 0.5 * dt) * relativity_factor / (pow(r, 3));*/
            pos[d] += vel[d] * 0.5 * dt;
        }
        pars[p].Par_SetPos(pos);
    } // for ( int p = 0; p < N_PARS; p++ )


    delete[] vel;
    delete[] acc;

    for (int d = 0; d < N_DIMS; d++)    delete[] force[d];
    delete[] force;

} // FUNCTION : Evolve_DKD



//--------------------------------------------------------------------------------
// Function    : Evolve_RK4
// Description : 
// Note        :
// Input       : 
// Output      : 
//--------------------------------------------------------------------------------
void Evolve_RK4( matrix &dens, particle *pars, const double dt )
{
} // FUNCTION : Evolve_RK4
