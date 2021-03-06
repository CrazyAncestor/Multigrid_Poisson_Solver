#include "classes.h"


matrix::matrix( int n, double hi )
{
    dim = n;
    h = hi;
    #if ( N_DIMS == 2 )
    cells = n*n;
    value = new double[n*n];
    #elif ( N_DIMS == 3)
    cells = n*n*n;
    value = new double[n*n*n];
    #endif

    // cell index increments
    did_x[0] = 1;       // i
    did_x[1] = dim;     // j
    did_x[2] = dim*dim; // k


    for ( int idx = 0; idx < cells; idx++ )    value[idx] = 0.0;

} // CONSTRUCTER : matrix::matrix( int n, double hi )



matrix::~matrix()
{
} // DESTRUCTURE : matrix::~matrix()



double absolute( double x )
{
    if ( x >= 0 )   return x;
    return -x;
} // FUNCTION : absolute()



void matrix::display()
{
    for ( int idx = 0; idx < cells; idx++ )
    {
        const int i = idx%dim;
        const int j = ( idx%(dim*dim) ) / dim;
        const int k = idx/(dim*dim);

        cout << value[idx] << " ";

        //if ( k == dim-1 )    cout << endl;
        if ( i == dim-1 )    cout << endl;

    } // for ( int idx = 0; idx < cells; idx++ )

} // FUNCTION : matrix::display()



void matrix::Error( const matrix &b )
{
    double sum = 0;
    double ave = 0;
    
    for ( int idx = 0; idx < cells; idx++ )
    {
        sum += absolute( this->value[idx] - b.value[idx] );
        ave += absolute( b.value[idx] );

    } // for ( int idx = 0; idx < cells; idx++ )

    cout << sum/ave << endl;

} //FUNCTION : matrix::Error

#   ifdef GPU
__global__
void SOR_smoothing_GPU(const bool even, const int dim, const double h, double* d_value, double* d_rho){

#   if(N_DIMS==2)
    const double frac = 0.25;
    const int i = blockIdx.x;
    const int j = threadIdx.x;
    int a;
    if (even) { a = 0; }
    else a = 1;
    if (i == 0 || i == dim - 1 || j == 0 || j == dim - 1 || (i + j) % 2 == a) {
    }
    else {
        int idx = i * dim + j;
        d_value[idx] += SOR_OMEGA * frac * (d_value[idx + 1] + d_value[idx - 1] +
            d_value[idx + dim] + d_value[idx - dim] -
            d_value[idx] * 4 - h * h * d_rho[idx]);
    }
#   elif(N_DIMS==3)
    const double frac = 1/6.;
    const int i = blockIdx.x;
    const int j = threadIdx.x / BOX_N;
    const int k = threadIdx.x % BOX_N;
    int a;
    if (even) { a = 0; }
    else a = 1;
    if (i == 0 || i == dim - 1 || j == 0 || j == dim - 1 || k == 0 || k == dim - 1 || (i + j + k) % 2 == a) {
    }
    else {
        int idx = i * dim * dim + j * dim + k;
        d_value[idx] += SOR_OMEGA * frac * (d_value[idx + 1] + d_value[idx - 1] +
            d_value[idx + dim] + d_value[idx - dim] + d_value[idx + dim * dim] + d_value[idx - dim * dim] -
            d_value[idx] * 6 - h * h * d_rho[idx]); 
    }
#   endif //#   if(N_DIMS==2)
} // FUNCTION :  SOR_smoothing_GPU
#endif

void matrix::SOR_smoothing(const matrix& rho, int steps)
{
#ifndef GPU
#   if ( N_DIMS == 2 )
    double fracs = 0.25;
    for (int t = 0; t < steps; t++) {
#   ifdef OMP_PARALLEL
#   pragma omp parallel num_threads(OMP_THREAD_NUM) 
        {
            const int tid = omp_get_thread_num();
#       pragma omp for collapse(1)
#   endif
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    if (i == 0 || i == dim - 1 || j == 0 || j == dim - 1 || (i + j) % 2 == 0) {
                    }
                    else {
                        int idx = i * dim + j;
                        this->value[idx] += SOR_OMEGA * fracs * (this->value[idx + 1] + this->value[idx - 1] +
                            this->value[idx + dim] + this->value[idx - dim] -
                            this->value[idx] * 4 - h * h * rho.value[idx]);
                    }

                }
            }
#   ifdef OMP_PARALLEL
#       pragma omp barrier
#       pragma omp for collapse(1)
#   endif
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    if (i == 0 || i == dim - 1 || j == 0 || j == dim - 1 || (i + j) % 2 == 1) {
                    }
                    else {
                        int idx = i * dim + j;
                        this->value[idx] += SOR_OMEGA * fracs * (this->value[idx + 1] + this->value[idx - 1] +
                            this->value[idx + dim] + this->value[idx - dim] -
                            this->value[idx] * 4 - h * h * rho.value[idx]);
                    }

                }
            }
#   ifdef OMP_PARALLEL
        }//pragma omp parallel num_threads(OMP_THREAD_NUM) 
#   endif//ifdef OMP_PARALLEL
        if (t % 1000 == 0)
        {
            //if (this->Residual_error(rho) < SOR_ERROR)break;
        }
    }//for (int t = 0; t < steps; t++)
#   elif ( N_DIMS == 3 )
    double fracs = 1 / 6.;
    for (int t = 0; t < steps; t++) {
#   ifdef OMP_PARALLEL
#   pragma omp parallel num_threads(OMP_THREAD_NUM) 
        {
            const int tid = omp_get_thread_num();
#       pragma omp for collapse(1)
#   endif
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    for (int k = 0; k < dim; k++) {
                        if (i == 0 || i == dim - 1 || j == 0 || j == dim - 1 || k == 0 || k == dim - 1 || (i + j + k) % 2 == 0) {
                        }
                        else {
                            int idx = i * dim * dim + j * dim + k;
                            this->value[idx] += SOR_OMEGA * fracs * (this->value[idx + 1] + this->value[idx - 1] +
                                this->value[idx + dim] + this->value[idx - dim] + this->value[idx + dim * dim] + this->value[idx - dim * dim] -
                                this->value[idx] * 6 - h * h * rho.value[idx]);
                        }
                    }

                }
            }
#   ifdef OMP_PARALLEL
#       pragma omp barrier
#       pragma omp for collapse(1)
#   endif
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    for (int k = 0; k < dim; k++) {
                        if (i == 0 || i == dim - 1 || j == 0 || j == dim - 1 || k == 0 || k == dim - 1 || (i + j + k) % 2 == 1) {
                        }
                        else {
                            int idx = i * dim * dim + j * dim + k;
                            this->value[idx] += SOR_OMEGA * fracs * (this->value[idx + 1] + this->value[idx - 1] +
                                this->value[idx + dim] + this->value[idx - dim] + this->value[idx + dim * dim] + this->value[idx - dim * dim] -
                                this->value[idx] * 6 - h * h * rho.value[idx]);
                        }
                    }

                }
            }
#   ifdef OMP_PARALLEL
        }//pragma omp parallel num_threads(OMP_THREAD_NUM) 
#   endif//ifdef OMP_PARALLEL
        if (t % 1000 == 0) 
        {
            //if (this->Residual_error(rho) < SOR_ERROR)break;
        }
}//for (int t = 0; t < steps; t++)
#   endif//if ( N_DIMS == 2 )

#   else//ifndef GPU
    double* d_value, * d_rho;
    int block_num = this->dim;
#   if ( N_DIMS == 2 )
    int block_scale = this->dim;
    int cell_num = pow(this->dim, 2);
#   elif ( N_DIMS == 3 )//if ( N_DIMS == 2 )
    int block_scale = pow(this->dim, 2);
    int cell_num = pow(this->dim, 3);
#   endif//if ( N_DIMS == 2 )
    cudaMalloc(&d_value, cell_num * sizeof(double));
    cudaMalloc(&d_rho, cell_num * sizeof(double));


    cudaMemcpy(d_value, this->value, cell_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_rho, rho.value, cell_num * sizeof(double), cudaMemcpyHostToDevice);

    for (int t = 0; t < steps; t++)
    {
        SOR_smoothing_GPU << < block_num, block_scale >> > (true, dim, h, d_value, d_rho);
        SOR_smoothing_GPU << < block_num, block_scale >> > (false, dim, h, d_value, d_rho);

        //if (t ==steps-1)
        {
            
        }
        
        //if (t % 1000 == 0)
        {
            /*cudaMemcpy(this->value, d_value, cells * sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(rho.value, d_rho, cells * sizeof(double), cudaMemcpyDeviceToHost);
            if (this->Residual_error(rho) < SOR_ERROR)break;*/
        }
    }
    cudaMemcpy(this->value, d_value, cell_num * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(rho.value, d_rho, cell_num * sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(d_value);//important
    cudaFree(d_rho);

#   endif//ifndef GPU
} // FUNCTION : matrix::SOR_smoothing




void matrix::SOR_Exact( const matrix &rho, int steps )
{
#ifndef GPU
#   if ( N_DIMS == 2 )
    double fracs = 0.25;
    for (int t = 0; t < steps; t++) {
#   ifdef OMP_PARALLEL
#   pragma omp parallel num_threads(OMP_THREAD_NUM) 
        {
            const int tid = omp_get_thread_num();
#       pragma omp for collapse(1)
#   endif
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    if (i == 0 || i == dim - 1 || j == 0 || j == dim - 1 || (i + j) % 2 == 0) {
                    }
                    else {
                        int idx = i * dim + j;
                        this->value[idx] += SOR_OMEGA * fracs * (this->value[idx + 1] + this->value[idx - 1] +
                            this->value[idx + dim] + this->value[idx - dim] -
                            this->value[idx] * 4 - h * h * rho.value[idx]);
                    }

                }
            }
#   ifdef OMP_PARALLEL
#       pragma omp barrier
#       pragma omp for collapse(1)
#   endif
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    if (i == 0 || i == dim - 1 || j == 0 || j == dim - 1 || (i + j) % 2 == 1) {
                    }
                    else {
                        int idx = i * dim + j;
                        this->value[idx] += SOR_OMEGA * fracs * (this->value[idx + 1] + this->value[idx - 1] +
                            this->value[idx + dim] + this->value[idx - dim] -
                            this->value[idx] * 4 - h * h * rho.value[idx]);
                    }

                }
            }
#   ifdef OMP_PARALLEL
            }//pragma omp parallel num_threads(OMP_THREAD_NUM) 
#   endif//ifdef OMP_PARALLEL
        if (t % 1000 == 0)
        {
            //if (this->Residual_error(rho) < SOR_ERROR)break;
        }
    }//for (int t = 0; t < steps; t++)
#   elif ( N_DIMS == 3 )
    double fracs = 1 / 6.;
    for (int t = 0; t < steps; t++) {
#   ifdef OMP_PARALLEL
#   pragma omp parallel num_threads(OMP_THREAD_NUM) 
        {
            const int tid = omp_get_thread_num();
#       pragma omp for collapse(1)
#   endif
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    for (int k = 0; k < dim; k++) {
                        if (i == 0 || i == dim - 1 || j == 0 || j == dim - 1 || k == 0 || k == dim - 1 || (i + j + k) % 2 == 0) {
                        }
                        else {
                            int idx = i * dim * dim + j * dim + k;
                            this->value[idx] += SOR_OMEGA * fracs * (this->value[idx + 1] + this->value[idx - 1] +
                                this->value[idx + dim] + this->value[idx - dim] + this->value[idx + dim * dim] + this->value[idx - dim * dim] -
                                this->value[idx] * 6 - h * h * rho.value[idx]);
}
                    }

                }
            }
#   ifdef OMP_PARALLEL
#       pragma omp barrier
#       pragma omp for collapse(1)
#   endif
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    for (int k = 0; k < dim; k++) {
                        if (i == 0 || i == dim - 1 || j == 0 || j == dim - 1 || k == 0 || k == dim - 1 || (i + j + k) % 2 == 1) {
                        }
                        else {
                            int idx = i * dim * dim + j * dim + k;
                            this->value[idx] += SOR_OMEGA * fracs * (this->value[idx + 1] + this->value[idx - 1] +
                                this->value[idx + dim] + this->value[idx - dim] + this->value[idx + dim * dim] + this->value[idx - dim * dim] -
                                this->value[idx] * 6 - h * h * rho.value[idx]);
                        }
                    }

                }
            }
#   ifdef OMP_PARALLEL
    }//pragma omp parallel num_threads(OMP_THREAD_NUM) 
#   endif//ifdef OMP_PARALLEL
        if (t % 1000 == 0)
        {
            //if (this->Residual_error(rho) < SOR_ERROR)break;
        }
}//for (int t = 0; t < steps; t++)
#   endif//if ( N_DIMS == 2 )

#   else//ifndef GPU
    double* d_value, * d_rho;
    int block_num = this->dim;
    #   if ( N_DIMS == 2 )
    int block_scale = this->dim;
    int cell_num = pow(this->dim, 2);
    #   elif ( N_DIMS == 3 )//if ( N_DIMS == 2 )
    int block_scale = pow(this->dim, 2);
    int cell_num = pow(this->dim, 3);
    #   endif//if ( N_DIMS == 2 )
    cudaMalloc(&d_value, cell_num * sizeof(double));
    cudaMalloc(&d_rho, cell_num * sizeof(double));


    cudaMemcpy(d_value, this->value, cell_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_rho, rho.value, cell_num * sizeof(double), cudaMemcpyHostToDevice);

    for (int t = 0; t < steps; t++)
    {
        SOR_smoothing_GPU << < block_num, block_scale >> > (true, dim, h, d_value, d_rho);
        SOR_smoothing_GPU << < block_num, block_scale >> > (false, dim, h, d_value, d_rho);

        //if (t ==steps-1)
        {

        }

        //if (t % 1000 == 0)
        {
            /*cudaMemcpy(this->value, d_value, cells * sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(rho.value, d_rho, cells * sizeof(double), cudaMemcpyDeviceToHost);
            if (this->Residual_error(rho) < SOR_ERROR)break;*/
        }
    }
    cudaMemcpy(this->value, d_value, cell_num * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(rho.value, d_rho, cell_num * sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(d_value);//important
    cudaFree(d_rho);

#   endif//ifndef GPU
} // FUNCTION : matrix::SOR_Exact


double matrix::averaging( int idx )
{
    // weight
    // 1 2 1 | 2 4 2 | 1 2 1
    // 2 4 2 | 4 8 4 | 2 4 2
    // 1 2 1 | 2 4 2 | 1 2 1
    int weight[3][3][3] = { { {1, 2, 1}, 
                              {2, 4, 2},
                              {1, 2, 1} }, 

                            { {2, 4, 2},
                              {4, 8, 4},
                              {2, 4, 2} },

                            { {1, 2, 1},
                              {2, 4, 2},
                              {1, 2, 1} } };

    const int i = idx%dim;
    const int j = ( idx%(dim*dim) ) / dim;
    const int k = idx/(dim*dim);

    if ( i == 0 )    for ( int dj = 0; dj < 3; dj++ )    for ( int dk = 0; dk < 3; dk++ )    weight[ 0][dj][dk] = 0;
    if ( j == 0 )    for ( int di = 0; di < 3; di++ )    for ( int dk = 0; dk < 3; dk++ )    weight[di][ 0][dk] = 0;
    if ( k == 0 )    for ( int di = 0; di < 3; di++ )    for ( int dj = 0; dj < 3; dj++ )    weight[di][dj][ 0] = 0;
    
    if ( dim%2 == 1 )
    {
        if ( i == dim-1 )    for ( int dj = 0; dj < 3; dj++ )    for ( int dk = 0; dk < 3; dk++ )    weight[ 2][dj][dk] = 0;
        if ( j == dim-1 )    for ( int di = 0; di < 3; di++ )    for ( int dk = 0; dk < 3; dk++ )    weight[di][ 2][dk] = 0;
        if ( k == dim-1 )    for ( int di = 0; di < 3; di++ )    for ( int dj = 0; dj < 3; dj++ )    weight[di][dj][ 2] = 0;
    } else // if ( dim%2 == 1 )
    {
        if ( i == dim-2 )    for ( int dj = 0; dj < 3; dj++ )    for ( int dk = 0; dk < 3; dk++ )   weight[ 2][dj][dk] = 0;
        if ( j == dim-2 )    for ( int di = 0; di < 3; di++ )    for ( int dk = 0; dk < 3; dk++ )   weight[di][ 2][dk] = 0;
        if ( k == dim-2 )    for ( int di = 0; di < 3; di++ )    for ( int dj = 0; dj < 3; dj++ )   weight[di][dj][ 2] = 0;
    } // if ( dim%2 == 1 ) ... else ...

    
    double avg = 0.0;
    double sum_w = 0.0;

    #if ( N_DIMS == 2 )
    for ( int di = -1; di < 2; di++ )
    {
        for ( int dj = -1; dj < 2; dj++ )
        {
            if ( weight[di+1][dj+1][1] == 0 )    continue;
            avg += weight[di+1][dj+1][1] * value[ idx + di*did_x[0] + dj*did_x[1] ]; 
            sum_w += weight[di+1][dj+1][1];
        }
    }
    #elif ( N_DIMS == 3 )
    for ( int di = -1; di < 2; di++ )
    {
        for ( int dj = -1; dj < 2; dj++ )
        {
            for ( int dk = -1; dk < 2; dk++ )
            {
                if ( weight[di+1][dj+1][dk+1] == 0 )    continue;
                avg += weight[di+1][dj+1][dk+1] * value[ idx + di*did_x[0] + dj*did_x[1] + dk*did_x[2] ]; 
                sum_w += weight[di+1][dj+1][dk+1];
            }
        }
    }
    #endif

    return avg / sum_w;

} // FUNCTION : matrix::averaging



matrix matrix::Restriction(){
    int dim2;
    if  ( dim%2 == 0 )  dim2 = dim / 2;
    else                dim2 = (dim+1) / 2;
    
    matrix r2h( dim2, h*2.0 );//important

    for ( int idx = 0; idx < cells; idx++ )
    {
        const int i = idx%dim;
        const int j = ( idx%(dim*dim) ) / dim;
        const int k = idx/(dim*dim);

        if( i%2 != 0 || j%2 != 0 || k%2 != 0 )    continue;
       
        const int idx2 = i/2 + dim2 * j/2 + dim2 * dim2 * k/2;
        r2h.value[idx2] = this->averaging( idx );

    } // for ( int idx = 0; idx < cells; idx++ )

    return r2h;

} // FUNCTION : matrix::Restriction



double matrix::insertion( int idx, int dim_in )
{
    int dim_new = dim_in;
    const int i = idx%dim_new;
    const int j = ( idx%(dim_new*dim_new) ) / dim_new;
    const int k = idx/(dim_new*dim_new);
    
    #if ( N_DIMS == 2 )
    if ( i == dim_new-1 && j == dim_new-1 )
    {
        if     ( i%2 == 0 && j%2 == 0 ) return   value[     i/2 + dim*    j/2 ];
        else if( i%2 == 1 && j%2 == 1 ) return   value[ (i-1)/2 + dim*(j-1)/2 ];
    }
    else if ( i == dim_new-1 )
    {
        if     ( i%2 == 0 && j%2 == 0 ) return   value[     i/2 + dim*    j/2 ];
        else if( i%2 == 0 && j%2 == 1 ) return ( value[     i/2 + dim*(j+1)/2 ] +
                                                 value[     i/2 + dim*(j-1)/2 ] ) * 0.5;
        else if( i%2 == 1 && j%2 == 0 ) return   value[ (i-1)/2 + dim*    j/2 ];
        else if( i%2 == 1 && j%2 == 1 ) return ( value[ (i-1)/2 + dim*(j+1)/2 ] +
                                                 value[ (i-1)/2 + dim*(j-1)/2 ] ) * 0.5;
    }
    else if ( j == dim_new-1 )
    {
        if     ( i%2 == 0 && j%2 == 0 ) return   value[     i/2 + dim*    j/2 ];
        else if( i%2 == 1 && j%2 == 0 ) return ( value[ (i+1)/2 + dim*    j/2 ] +
                                                 value[ (i-1)/2 + dim*    j/2 ] ) * 0.5;
        else if( i%2 == 0 && j%2 == 1 ) return   value[     i/2 + dim*(j-1)/2 ];
        else if( i%2 == 1 && j%2 == 1 ) return ( value[ (i+1)/2 + dim*(j-1)/2 ] +
                                                 value[ (i-1)/2 + dim*(j-1)/2 ] ) * 0.5;
    }
    else
    {
        if     ( i%2 == 0 && j%2 == 0 ) return   value[     i/2 + dim*    j/2 ];
        else if( i%2 == 0 && j%2 == 1 ) return ( value[     i/2 + dim*(j+1)/2 ] +
                                                 value[     i/2 + dim*(j-1)/2 ] ) * 0.5;
        else if( i%2 == 1 && j%2 == 0 ) return ( value[ (i+1)/2 + dim*    j/2 ] +
                                                 value[ (i-1)/2 + dim*    j/2 ] ) * 0.5;
        else if( i%2 == 1 && j%2 == 1 ) return ( value[ (i+1)/2 + dim*(j+1)/2 ] +
                                                 value[ (i+1)/2 + dim*(j-1)/2 ] +
                                                 value[ (i-1)/2 + dim*(j+1)/2 ] +
                                                 value[ (i-1)/2 + dim*(j-1)/2 ] ) * 0.25;
    }

    #elif ( N_DIMS == 3 )
    if ( i == dim_new-1 && j == dim_new-1 && k == dim_new-1 )
    {
        if      ( i%2 == 0 && j%2 == 0 && k%2 == 0 )    return   value[     i/2 + dim*    j/2 + dim*dim*    k/2 ];
        else if ( i%2 == 1 && j%2 == 1 && k%2 == 1 )    return   value[ (i-1)/2 + dim*(j-1)/2 + dim*dim*(k-1)/2 ];
    }
    else if ( i == dim_new-1 && j == dim_new-1 )
    {
        if      ( i%2 == 0 && j%2 == 0 && k%2 == 0 )    return   value[     i/2 + dim*    j/2 + dim*dim*    k/2 ];
        else if ( i%2 == 0 && j%2 == 0 && k%2 == 1 )    return ( value[     i/2 + dim*    j/2 + dim*dim*(k+1)/2 ] + 
                                                                 value[     i/2 + dim*    j/2 + dim*dim*(k-1)/2 ] ) * 0.5;
        else if ( i%2 == 1 && j%2 == 1 && k%2 == 0 )    return   value[ (i-1)/2 + dim*(j-1)/2 + dim*dim*    k/2 ];
        else if ( i%2 == 1 && j%2 == 1 && k%2 == 1 )    return ( value[ (i-1)/2 + dim*(j-1)/2 + dim*dim*(k+1)/2 ] + 
                                                                 value[ (i-1)/2 + dim*(j-1)/2 + dim*dim*(k-1)/2 ] ) * 0.5;
    }
    else if ( i == dim_new-1 && k == dim_new-1 )
    {
        if      ( i%2 == 0 && j%2 == 0 && k%2 == 0 )    return   value[     i/2 + dim*    j/2 + dim*dim*    k/2 ];
        else if ( i%2 == 0 && j%2 == 1 && k%2 == 0 )    return ( value[     i/2 + dim*(j+1)/2 + dim*dim*    k/2 ] +
                                                                 value[     i/2 + dim*(j-1)/2 + dim*dim*    k/2 ] ) * 0.5;
        else if ( i%2 == 1 && j%2 == 0 && k%2 == 1 )    return   value[ (i-1)/2 + dim*    j/2 + dim*dim*(k-1)/2 ];
        else if ( i%2 == 1 && j%2 == 1 && k%2 == 1 )    return ( value[ (i-1)/2 + dim*(j+1)/2 + dim*dim*(k-1)/2 ] +
                                                                 value[ (i-1)/2 + dim*(j-1)/2 + dim*dim*(k-1)/2 ] ) * 0.5;
    }
    else if ( j == dim_new-1 && k == dim_new-1 )
    {
        if      ( i%2 == 0 && j%2 == 0 && k%2 == 0 )    return   value[     i/2 + dim*    j/2 + dim*dim*    k/2 ];
        else if ( i%2 == 0 && j%2 == 1 && k%2 == 1 )    return   value[     i/2 + dim*(j-1)/2 + dim*dim*(k-1)/2 ];
        else if ( i%2 == 1 && j%2 == 0 && k%2 == 0 )    return ( value[ (i+1)/2 + dim*    j/2 + dim*dim*    k/2 ] +
                                                                 value[ (i-1)/2 + dim*    j/2 + dim*dim*    k/2 ] ) * 0.5;
        else if ( i%2 == 1 && j%2 == 1 && k%2 == 1 )    return ( value[ (i+1)/2 + dim*(j-1)/2 + dim*dim*(k-1)/2 ] +
                                                                 value[ (i-1)/2 + dim*(j-1)/2 + dim*dim*(k-1)/2 ] ) * 0.5;
    }
    else if ( i == dim_new-1 )
    {
        if      ( i%2 == 0 && j%2 == 0 && k%2 == 0 )    return   value[     i/2 + dim*    j/2 + dim*dim*    k/2 ];
        else if ( i%2 == 0 && j%2 == 0 && k%2 == 1 )    return ( value[     i/2 + dim*    j/2 + dim*dim*(k+1)/2 ] +
                                                                 value[     i/2 + dim*    j/2 + dim*dim*(k-1)/2 ] ) * 0.5;
        else if ( i%2 == 0 && j%2 == 1 && k%2 == 0 )    return ( value[     i/2 + dim*(j+1)/2 + dim*dim*    k/2 ] + 
                                                                 value[     i/2 + dim*(j-1)/2 + dim*dim*    k/2 ] ) * 0.5;
        else if ( i%2 == 0 && j%2 == 1 && k%2 == 1 )    return ( value[     i/2 + dim*(j+1)/2 + dim*dim*(k+1)/2 ] +
                                                                 value[     i/2 + dim*(j+1)/2 + dim*dim*(k-1)/2 ] +
                                                                 value[     i/2 + dim*(j-1)/2 + dim*dim*(k+1)/2 ] +
                                                                 value[     i/2 + dim*(j-1)/2 + dim*dim*(k-1)/2 ] ) * 0.25;
        else if ( i%2 == 1 && j%2 == 0 && k%2 == 0 )    return   value[ (i-1)/2 + dim*    j/2 + dim*dim*    k/2 ];
        else if ( i%2 == 1 && j%2 == 0 && k%2 == 1 )    return ( value[ (i-1)/2 + dim*    j/2 + dim*dim*(k+1)/2 ] +
                                                                 value[ (i-1)/2 + dim*    j/2 + dim*dim*(k-1)/2 ] ) * 0.5;
        else if ( i%2 == 1 && j%2 == 1 && k%2 == 0 )    return ( value[ (i-1)/2 + dim*(j+1)/2 + dim*dim*    k/2 ] +
                                                                 value[ (i-1)/2 + dim*(j-1)/2 + dim*dim*    k/2 ] ) * 0.5;
        else if ( i%2 == 1 && j%2 == 1 && k%2 == 1 )    return ( value[ (i-1)/2 + dim*(j+1)/2 + dim*dim*(k+1)/2 ] +
                                                                 value[ (i-1)/2 + dim*(j+1)/2 + dim*dim*(k-1)/2 ] +
                                                                 value[ (i-1)/2 + dim*(j-1)/2 + dim*dim*(k+1)/2 ] +
                                                                 value[ (i-1)/2 + dim*(j-1)/2 + dim*dim*(k-1)/2 ] ) * 0.25;
    }
    else if ( j == dim_new-1 )
    {
        if      ( i%2 == 0 && j%2 == 0 && k%2 == 0 )    return   value[     i/2 + dim*    j/2 + dim*dim*    k/2 ];
        else if ( i%2 == 0 && j%2 == 0 && k%2 == 1 )    return ( value[     i/2 + dim*    j/2 + dim*dim*(k+1)/2 ] +
                                                                 value[     i/2 + dim*    j/2 + dim*dim*(k-1)/2 ] ) * 0.5;
        else if ( i%2 == 0 && j%2 == 1 && k%2 == 0 )    return   value[     i/2 + dim*(j-1)/2 + dim*dim*    k/2 ]; 
        else if ( i%2 == 0 && j%2 == 1 && k%2 == 1 )    return ( value[     i/2 + dim*(j-1)/2 + dim*dim*(k+1)/2 ] + 
                                                                 value[     i/2 + dim*(j-1)/2 + dim*dim*(k-1)/2 ] ) * 0.5;
        else if ( i%2 == 1 && j%2 == 0 && k%2 == 0 )    return ( value[ (i+1)/2 + dim*    j/2 + dim*dim*(k+1)/2 ] + 
                                                                 value[ (i-1)/2 + dim*    j/2 + dim*dim*(k-1)/2 ] ) * 0.5;
        else if ( i%2 == 1 && j%2 == 0 && k%2 == 1 )    return ( value[ (i+1)/2 + dim*    j/2 + dim*dim*(k+1)/2 ] +
                                                                 value[ (i+1)/2 + dim*    j/2 + dim*dim*(k-1)/2 ] +
                                                                 value[ (i-1)/2 + dim*    j/2 + dim*dim*(k+1)/2 ] +
                                                                 value[ (i-1)/2 + dim*    j/2 + dim*dim*(k-1)/2 ] ) * 0.25;
        else if ( i%2 == 1 && j%2 == 1 && k%2 == 0 )    return ( value[ (i+1)/2 + dim*(j-1)/2 + dim*dim*    k/2 ] +
                                                                 value[ (i-1)/2 + dim*(j-1)/2 + dim*dim*    k/2 ] ) * 0.5;
        else if ( i%2 == 1 && j%2 == 1 && k%2 == 1 )    return ( value[ (i+1)/2 + dim*(j-1)/2 + dim*dim*(k+1)/2 ] +
                                                                 value[ (i+1)/2 + dim*(j-1)/2 + dim*dim*(k-1)/2 ] +
                                                                 value[ (i-1)/2 + dim*(j-1)/2 + dim*dim*(k+1)/2 ] +
                                                                 value[ (i-1)/2 + dim*(j-1)/2 + dim*dim*(k-1)/2 ] ) * 0.25;
    }
    else if ( k == dim_new-1 )
    {
        if      ( i%2 == 0 && j%2 == 0 && k%2 == 0 )    return   value[     i/2 + dim*    j/2 + dim*dim*    k/2 ];
        else if ( i%2 == 0 && j%2 == 0 && k%2 == 1 )    return   value[     i/2 + dim*    j/2 + dim*dim*(k-1)/2 ];
        else if ( i%2 == 0 && j%2 == 1 && k%2 == 0 )    return ( value[     i/2 + dim*(j+1)/2 + dim*dim*    k/2 ] +
                                                                 value[     i/2 + dim*(j-1)/2 + dim*dim*    k/2 ] ) * 0.5;
        else if ( i%2 == 0 && j%2 == 1 && k%2 == 1 )    return ( value[     i/2 + dim*(j+1)/2 + dim*dim*(k-1)/2 ] + 
                                                                 value[     i/2 + dim*(j-1)/2 + dim*dim*(k-1)/2 ] ) * 0.5;
        else if ( i%2 == 1 && j%2 == 0 && k%2 == 0 )    return ( value[ (i+1)/2 + dim*    j/2 + dim*dim*    k/2 ] + 
                                                                 value[ (i-1)/2 + dim*    j/2 + dim*dim*    k/2 ] ) * 0.5;
        else if ( i%2 == 1 && j%2 == 0 && k%2 == 1 )    return ( value[ (i+1)/2 + dim*    j/2 + dim*dim*    k/2 ] + 
                                                                 value[ (i-1)/2 + dim*    j/2 + dim*dim*(k-1)/2 ] ) * 0.5;
        else if ( i%2 == 1 && j%2 == 1 && k%2 == 0 )    return ( value[ (i+1)/2 + dim*(j+1)/2 + dim*dim*    k/2 ] +
                                                                 value[ (i+1)/2 + dim*(j-1)/2 + dim*dim*    k/2 ] +
                                                                 value[ (i-1)/2 + dim*(j+1)/2 + dim*dim*    k/2 ] +
                                                                 value[ (i-1)/2 + dim*(j-1)/2 + dim*dim*    k/2 ] ) * 0.25;
        else if ( i%2 == 1 && j%2 == 1 && k%2 == 1 )    return ( value[ (i+1)/2 + dim*(j+1)/2 + dim*dim*(k-1)/2 ] + 
                                                                 value[ (i+1)/2 + dim*(j-1)/2 + dim*dim*(k-1)/2 ] +
                                                                 value[ (i-1)/2 + dim*(j+1)/2 + dim*dim*(k-1)/2 ] +
                                                                 value[ (i-1)/2 + dim*(j-1)/2 + dim*dim*(k-1)/2 ] ) * 0.25;
    }
    else 
    {
        if      ( i%2 == 0 && j%2 == 0 && k%2 == 0 )    return   value[     i/2 + dim*    j/2 + dim*dim*    k/2 ];
        else if ( i%2 == 0 && j%2 == 0 && k%2 == 1 )    return ( value[     i/2 + dim*    j/2 + dim*dim*(k+1)/2 ] +
                                                                 value[     i/2 + dim*    j/2 + dim*dim*(k-1)/2 ] ) * 0.5;
        else if ( i%2 == 0 && j%2 == 1 && k%2 == 0 )    return ( value[     i/2 + dim*(j+1)/2 + dim*dim*    k/2 ] +
                                                                 value[     i/2 + dim*(j-1)/2 + dim*dim*    k/2 ] ) * 0.5;
        else if ( i%2 == 0 && j%2 == 1 && k%2 == 1 )    return ( value[     i/2 + dim*(j+1)/2 + dim*dim*(k+1)/2 ] +
                                                                 value[     i/2 + dim*(j+1)/2 + dim*dim*(k-1)/2 ] +
                                                                 value[     i/2 + dim*(j-1)/2 + dim*dim*(k+1)/2 ] +
                                                                 value[     i/2 + dim*(j-1)/2 + dim*dim*(k-1)/2 ] ) * 0.25;
        else if ( i%2 == 1 && j%2 == 0 && k%2 == 0 )    return ( value[ (i+1)/2 + dim*    j/2 + dim*dim*    k/2 ] +
                                                                 value[ (i-1)/2 + dim*    j/2 + dim*dim*    k/2 ] ) * 0.5;
        else if ( i%2 == 1 && j%2 == 0 && k%2 == 1 )    return ( value[ (i+1)/2 + dim*    j/2 + dim*dim*(k+1)/2 ] +
                                                                 value[ (i+1)/2 + dim*    j/2 + dim*dim*(k-1)/2 ] +
                                                                 value[ (i-1)/2 + dim*    j/2 + dim*dim*(k+1)/2 ] +
                                                                 value[ (i-1)/2 + dim*    j/2 + dim*dim*(k-1)/2 ] ) * 0.25;
        else if ( i%2 == 1 && j%2 == 1 && k%2 == 0 )    return ( value[ (i+1)/2 + dim*(j+1)/2 + dim*dim*    k/2 ] +
                                                                 value[ (i+1)/2 + dim*(j-1)/2 + dim*dim*    k/2 ] +
                                                                 value[ (i-1)/2 + dim*(j+1)/2 + dim*dim*    k/2 ] +
                                                                 value[ (i-1)/2 + dim*(j-1)/2 + dim*dim*    k/2 ] ) * 0.25;
        else if ( i%2 == 1 && j%2 == 1 && k%2 == 1 )    return ( value[ (i+1)/2 + dim*(j+1)/2 + dim*dim*(k+1)/2 ] +
                                                                 value[ (i+1)/2 + dim*(j+1)/2 + dim*dim*(k-1)/2 ] +
                                                                 value[ (i+1)/2 + dim*(j-1)/2 + dim*dim*(k+1)/2 ] +
                                                                 value[ (i+1)/2 + dim*(j-1)/2 + dim*dim*(k-1)/2 ] +
                                                                 value[ (i-1)/2 + dim*(j+1)/2 + dim*dim*(k+1)/2 ] +
                                                                 value[ (i-1)/2 + dim*(j+1)/2 + dim*dim*(k-1)/2 ] +
                                                                 value[ (i-1)/2 + dim*(j-1)/2 + dim*dim*(k+1)/2 ] +
                                                                 value[ (i-1)/2 + dim*(j-1)/2 + dim*dim*(k-1)/2 ] ) * 0.125;
    }
    #endif

} // FUNCTION : matrix::insertion



matrix matrix::Interpolation( int dim_in )
{
    int dim_new = dim_in;
    #if ( N_DIMS == 2 )
    int cells_new = dim_new * dim_new;
    #elif ( N_DIMS == 3 )
    int cells_new = dim_new * dim_new * dim_new;
    #endif
    matrix r ( dim_new, h/2.0 );
    
    for ( int idx = 0; idx < cells_new; idx++ )    r.value[idx] = this->insertion( idx, dim_new );
    
    return r;

} // FUNCTION : matrix::Interpolation

double matrix::Residual_error(const matrix& rho)
{
    double sum = 0.;
    for (int idx = 0; idx < cells; idx++)
    {
        const int i = idx % dim;
        const int j = (idx % (dim * dim)) / dim;
        const int k = idx / (dim * dim);

#if ( N_DIMS == 2 )
        if (i == 0 || i == dim - 1 || j == 0 || j == dim - 1)    continue;
        sum += 0.25 * (this->value[idx + did_x[0]] + this->value[idx - did_x[0]] +
            this->value[idx + did_x[1]] + this->value[idx - did_x[1]] -
            this->value[idx] * 4 - h * h * rho.value[idx]) / h / h / rho.value[idx];
#elif ( N_DIMS == 3 )
        if (i == 0 || i == dim - 1 || j == 0 || j == dim - 1 || k == 0 || k == dim - 1)    continue;
        sum += (this->value[idx + did_x[0]] + this->value[idx - did_x[0]] +
            this->value[idx + did_x[1]] + this->value[idx - did_x[1]] +
            this->value[idx + did_x[2]] + this->value[idx - did_x[2]] -
            this->value[idx] * 6 - h * h * rho.value[idx]) / h / h  / 6./ rho.value[idx];
#endif

    } // for( int idx = 0; idx < cells; idx++ ) 

    return sum;

} // FUNCTION : double matrix::Residual

matrix matrix::Residual( const matrix & rho )
{
    matrix res( this->dim, this->h );
    for (int idx = 0; idx < cells; idx++)
    {
        const int i = idx % dim;
        const int j = (idx % (dim * dim)) / dim;
        const int k = idx / (dim * dim);

#if ( N_DIMS == 2 )
        if (i == 0 || i == dim - 1 || j == 0 || j == dim - 1)  res.value[idx] = 0;//continue;//important-bc.problem.
        else {
            res.value[idx] = 0.25 * (this->value[idx + did_x[0]] + this->value[idx - did_x[0]] +
                this->value[idx + did_x[1]] + this->value[idx - did_x[1]] -
                this->value[idx] * 4 - h * h * rho.value[idx]) / h / h;
        }
#elif ( N_DIMS == 3 )
        if (i == 0 || i == dim - 1 || j == 0 || j == dim - 1 || k == 0 || k == dim - 1)  res.value[idx] = 0;  //continue;//important-bc.problem.
        else {
            res.value[idx] = (this->value[idx + did_x[0]] + this->value[idx - did_x[0]] +
                this->value[idx + did_x[1]] + this->value[idx - did_x[1]] +
                this->value[idx + did_x[2]] + this->value[idx - did_x[2]] -
                this->value[idx] * 6 - h * h * rho.value[idx] ) / h / h / 6.;
        }
        #endif
        
    } // for( int idx = 0; idx < cells; idx++ ) 

    return res;

} // FUNCTION : matrix matrix::Residual



matrix matrix::Laplacian()
{
    matrix lap( this->dim, this->h );
    
    for( int idx = 0; idx < cells; idx++ ) 
    {
        const int i = idx%dim;
        const int j = ( idx%(dim*dim) ) / dim;
        const int k = idx/(dim*dim);
        
        #if ( N_DIMS == 2 )
        if ( i == 0 || i == dim-1 || j == 0 || j == dim-1 )    continue;
        lap.value[idx] = 0.25 * ( this->value[ idx+did_x[0] ] + this->value[ idx-did_x[0] ] + 
                                  this->value[ idx+did_x[1] ] + this->value[ idx-did_x[1] ] - 
                                  this->value[ idx ]*4  ) / h / h;
        #elif ( N_DIMS == 3 )
        if ( i == 0 || i == dim-1 || j == 0 || j == dim-1 || k == 0 || k == dim-1 )    continue;
        lap.value[idx] = ( this->value[ idx+did_x[0] ] + this->value[ idx-did_x[0] ] + 
                           this->value[ idx+did_x[1] ] + this->value[ idx-did_x[1] ] + 
                           this->value[ idx+did_x[2] ] + this->value[ idx-did_x[2] ] - 
                           this->value[ idx ]*6  ) / h / h / 6.;
        #endif
        
    } // for( int idx = 0; idx < cells; idx++ ) 

    return lap;

} // FUNCTION : matrix::Laplacian



void matrix::init_density()
{
    for( int idx = 0; idx < cells; idx++ ) 
    {
        const int i = idx%dim;
        const int j = ( idx%(dim*dim) ) / dim;
        const int k = idx/(dim*dim);
        
        double x = h*i;
        double y = h*j;
        double z = h*k;
        
        this->value[idx] = -2.*sin(x)*sin(y)*sin(z);
    } // for( int idx = 0; idx < cells; idx++ ) 

} //FUNCTION : matrix::init_density



void matrix::init_potential()
{
    for( int idx = 0; idx < cells; idx++ )    this->value[idx] = BG_POTENTIAL;
} //FUNCTION : matrix::init_density



void matrix::reset()
{
    for( int idx = 0; idx < cells; idx++ )    this->value[idx] = 0.0;
    
} // FUNCTION : matrix::reset_particle_density



double matrix::get_h()
{
    return this->h;
} // FUNCTION : matrix::get_h



double matrix::get_dim()
{
    return this->dim;
} // FUNCTION : matrix::get_dim



double matrix::get_value( int idx )
{
    return this->value[idx];
} // FUNCTION : matrix::get_value



void matrix::set_value( int idx, double val )
{
    this->value[idx] = val;
} // FUNCTION : matrix::set_value



void matrix::add_value( int idx, double val )
{
    this->value[idx] += val;
} // FUNCTION : matrix::add_value



void matrix::input_answer( int idx, double ans )
{
    this->value[idx] = ans;
} // FUNCTION : matrix::input_answer
  


matrix matrix::operator+( const matrix &b )
{
    matrix tmp( dim, h );

    for ( int idx = 0; idx < cells; idx++ )    tmp.value[idx] = this->value[idx] + b.value[idx];

    return tmp;

} // FUNCTION : matrix::operator+



matrix matrix::operator-( const matrix &b )
{
    matrix tmp( dim, h );

    for ( int idx = 0; idx < cells; idx++ )    tmp.value[idx] = this->value[idx] - b.value[idx];
    
    return tmp;

} // FUNCTION : matrix::operator-

void  matrix::set_boundary(double bc) {
    for (int idx = 0; idx < cells; idx++)
    {
        const int i = idx % dim;
        const int j = (idx % (dim * dim)) / dim;
        const int k = idx / (dim * dim);

#if ( N_DIMS == 2 )
        if (i == 0 || i == dim - 1 || j == 0 || j == dim - 1)  this->value[idx] = bc;
        
#elif ( N_DIMS == 3 )
        if (i == 0 || i == dim - 1 || j == 0 || j == dim - 1 || k == 0 || k == dim - 1)  this->value[idx] = bc;
        
#endif

    } // for( int idx = 0; idx < cells; idx++ ) 
} // FUNCTION : matrix::set_boundary



