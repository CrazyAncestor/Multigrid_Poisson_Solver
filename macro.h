// constants 
#define PI              3.14159265359


// box
#define BOX_DX          BOX_L / (BOX_N-1)       // NOT cell center grid

#if ( N_DIMS == 2 )
#define N_CELLS         BOX_N*BOX_N
#elif ( N_DIMS == 3 )
#define N_CELLS         BOX_N*BOX_N*BOX_N
#endif


// problem number
#define PROB_SINWAVE    1
#define PROB_TWOBODY    2
#define PROB_NBODY      3
#define PROB_SOLARORBIT 4


// solver type
#define SOR             1
#define V_CYCLE         2
#define W_CYCLE         3
#define FAS             4
#define FMG             5


// particle mass to cell mass
#define NGP             1                       // Nearest-Grid-Point
#define CIC             2                       // Cloud-in-Cell
#define TSC             3                       // Triangular-Shape-Cloud


// particle moving scheme
#define EULER           1
#define KDK             2
#define DKD             3
#define RK4             4

