#include "classes.h"

#include "solver.h"


//--------------------------------------------------------------------------------
// Function    : Solver_Potential
// Description : 
// Note        :
// Input       : phi  : Potential of the test problem
//               dens : Density
// Output      : Solved potential 2D array
//--------------------------------------------------------------------------------
matrix Solver_Potential( matrix pot, matrix dens )
{
#   if ( POT_SOLVER == SOR ) 
    matrix solution = SOR_Method( pot, dens );
#   elif ( POT_SOLVER == V_CYCLE ) 
    matrix solution = V_Cycle( pot, dens );
#   elif ( POT_SOLVER == W_CYCLE ) 
    matrix solution = W_Cycle( pot, dens );
#   elif ( POT_SOLVER == FAS ) 
    // empty now
#   elif ( POT_SOLVER == FMG ) 
    matrix solution = FMG_Method( pot, dens, 1 );
#   endif // #if ( POT_SOLVER == ... )

    return solution;

} // FUNCTION : Solver_Potential



//--------------------------------------------------------------------------------
// Function    : Solver_Force
// Description : 
// Note        :
// Input       : pot  : Potential of the test problem
// Output      : Solved force 3D array
//--------------------------------------------------------------------------------
void Solver_Force( matrix pot, double **force )
{
    const double dx = BOX_DX;
    const double softening = 0.01;
    double r;
    int di, dj, dk;

    const int did_x[3] = { 1, BOX_N, BOX_N*BOX_N };
    for ( int d = 0; d < N_DIMS; d++ )
    {
        if ( d == 0 )
        {
            di = did_x[d];
            dj = 0;
            dk = 0;
        } else if ( d == 1 )
        {
            di = 0;
            dj = did_x[d];
            dk = 0;
        } else if ( d == 2 )
        {
            di = 0;
            dj = 0;
            dk = did_x[d];
        }// if ( d == 0 ) ... else if ...

        for ( int idx = 0; idx < N_CELLS; idx++ )
        {
            const int i = idx%BOX_N;
            const int j = (idx%(BOX_N*BOX_N)) / BOX_N;
            const int k = idx/(BOX_N*BOX_N);

            if ( (d == 0 && i == 0) || (d == 1 && j == 0) || (d == 2 && k == 0) )
            {
                force[d][idx] = -( pot.get_value( idx+di+dj+dk ) - pot.get_value( idx ) ) / dx;
            } else if ( (d == 0 && i == BOX_N-1) || ( d == 1 && j == BOX_N-1) || (d == 2 && k == BOX_N-1) )
            {
                force[d][idx] = -( pot.get_value( idx ) - pot.get_value( idx-di-dj-dk ) ) / dx;
            } else
            {
                r = 1 * pow((pow((i - 0.5 * BOX_N), 2) + pow((j - 0.5 * BOX_N), 2) +pow((k - 0.5 * BOX_N), 2) + pow(softening, 2)), 0.5);
                force[d][idx] = -0.5 * (pot.get_value(idx + di + dj + dk) - pot.get_value(idx - di - dj - dk)) / dx;
                //force[d][idx] = - 100 * pow(r, -2);
                //force[d][idx] = - 1 * pow(r,-2);
            } // if ( i == 0 || j == 0 ) ... else if ... else ...
            
        } // for ( int idx = 0; idx < N_CELLS; idx++ )

    } // for ( int d = 0; d < N_DIMS; d++ )

} // FUNCTION : Solver_Force



//--------------------------------------------------------------------------------
// Function    : SOR_Method
// Description : 
// Note        :
// Input       : phi  : Potential of the test problem
//               dens : Density
// Output      : Solved potential 2D array
//--------------------------------------------------------------------------------
matrix SOR_Method( matrix phi, matrix dens )
{   
    phi.SOR_Exact( dens, SOR_EXACT_STEP );
    
    return phi;

} // FUNCTION : V_Cycle



//--------------------------------------------------------------------------------
// Function    : V_Cycle
// Description : 
// Note        :
// Input       : phi  : Potential of the test problem
//               dens : Density
// Output      : Solved potential 2D array
//--------------------------------------------------------------------------------
matrix V_Cycle( matrix phi, matrix dens )
{   
    const int smooth_step  = SOR_SMOOTH_STEP;
    const int exact_step   = SOR_EXACT_STEP;

    const int n = dens.get_dim();

    // the smallest grid size, do the exact solver
    if ( n <= 3 )
    {
          // Exact solver
          phi.SOR_smoothing( dens, exact_step );
    }

    else{
        // Pre-Smoothing
        phi.SOR_smoothing( dens, smooth_step );
        
        // Compute Residual Errors
        matrix r = phi.Residual( dens );
    
        // Restriction
        matrix rhs = r.Restriction();
        
        matrix eps( rhs.get_dim(), rhs.get_h() );
        
        // Go to smaller grid size
        eps = V_Cycle( eps, rhs );  

        // Prolongation and Correction
        matrix phi_corr = eps.Interpolation(phi.get_dim());//important
        phi_corr.set_boundary(0.);//important
        phi = phi + phi_corr;//important
        
        // Post-Smoothing
    }
    
    phi.SOR_smoothing( dens, smooth_step );
    
    return phi;

} // FUNCTION : V_Cycle



//--------------------------------------------------------------------------------
// Function    : W_Cycle
// Description : 
// Note        :
// Input       : phi  : Potential of the test problem
//               dens : Density
//               LR   : 0 : the left valley of the W cycle
//                      1 : the right valley of the W cycle
// Output      : Solved potential 2D array
//--------------------------------------------------------------------------------
matrix W_Cycle( matrix phi, matrix dens )
{
     const int smooth_step  = SOR_SMOOTH_STEP;
     const int exact_step   = SOR_EXACT_STEP;
     const int n = dens.get_dim();
     
     // the smallest grid size, do the exact solver
     if ( n <= 3 )
     {
          // Exact solver
          phi.SOR_smoothing(dens, exact_step);
     }
     // stop recursion at smallest grid size
     
     else // if ( n <= 3 ) 
     {
          // Pre-Smoothing
          phi.SOR_smoothing( dens, smooth_step );
          
          // Restriction
          matrix r   = phi.Residual( dens );
          matrix rhs = r.Restriction();
          matrix eps( rhs.get_dim(), rhs.get_h() );
          
          // Left valley
          eps = W_Cycle( eps, rhs );
  
          // Right valley
          eps = W_Cycle( eps, rhs );
          
          // Prolongation
          matrix phi_corr = eps.Interpolation(phi.get_dim());//important
          phi_corr.set_boundary(0.);//important
          phi = phi + phi_corr;//important
          

          // Post-Smoothing
          phi.SOR_smoothing( dens, smooth_step );
          
     } // else //if ( n <= 3 ) 
     
     return phi;
   
} // FUNCTION : W_Cycle



//--------------------------------------------------------------------------------
// Function    : FAS_Method
// Description : 
// Note        :
// Input       : phi  : Potential of the test problem
//               dens : Density
// Output      : Solved potential 2D array
//--------------------------------------------------------------------------------
matrix FAS_Method( matrix phi, matrix dens )
{   
    const int smooth_step  = SOR_SMOOTH_STEP;
    const int exact_step   = SOR_EXACT_STEP;
    const int n = dens.get_dim();

    // the smallest grid size, do the exact solver
    if ( n <= 3 )
    {
          // Exact solver
          matrix eps( phi.get_dim(), phi.get_h() );
          matrix dens2 = dens.Restriction();
          eps.SOR_smoothing( dens2, exact_step );
    }

    else{
        // Pre-Smoothing
        phi.SOR_smoothing( dens, smooth_step );
        
        // Compute Residual Errors
        matrix res = phi.Residual( dens );
    
        // Restriction
        matrix res_coarse = res.Restriction();
        
        matrix phi_coarse = phi.Restriction();

        matrix dens_coarse = phi_coarse.Laplacian() +res_coarse;
        
        // Go to smaller grid size
        phi_coarse = V_Cycle( phi_coarse, dens_coarse );  
    
        // Prolongation and Correction
        matrix corr = phi_coarse -phi.Restriction();
        phi = phi + corr.Interpolation( phi.get_dim() );
        
        // Post-Smoothing
        phi.SOR_smoothing( dens, smooth_step );
    }
    
    return phi;

} // FUNCTION : FAS_Method



//--------------------------------------------------------------------------------
// Function    : FMG_Method
// Description : 
// Note        :
// Input       : phi  : Potential of the test problem
//               dens : Density
//               n_fmg: On each working level, one applies n_fmg MG cycles
// Output      : Solved potential 2D array
//--------------------------------------------------------------------------------
matrix FMG_Method( matrix phi, matrix dens, int n_fmg)
{   
    const int smooth_step  = SOR_SMOOTH_STEP;
    const int exact_step   = SOR_EXACT_STEP;
    const int n = dens.get_dim();

    matrix phi_coarse = phi.Restriction();
    matrix dens_coarse = dens.Restriction();
    
   /* //phi = phi_coarse.Interpolation(n);
    phi.display();
    phi = V_Cycle(phi, dens);
    phi.display();*/
    if(n<3){
        phi.SOR_smoothing( dens, exact_step );
    }
    else{
        matrix phi_coarse = phi.Restriction();
        matrix dens_coarse = dens.Restriction();
        phi_coarse = FMG_Method( phi_coarse, dens_coarse,n_fmg);
        phi = phi_coarse.Interpolation(n);
        phi.set_boundary(BG_POTENTIAL);//important-bc.problem
        for(int i=0;i<n_fmg;i++){
            phi = V_Cycle( phi, dens );
        }
    }
    return phi;

} // FUNCTION : FMG_Method
