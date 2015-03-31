extern int g_Ngrid_max;
//extern double sigma;
extern int g_is_ball;
extern int g_n_modes;
extern const Target_1dim* g_targ_1d;
extern const Ndim_histogram* g_targp;
extern Target_peak_1dim* g_peaks;
extern const gsl_rng_type* g_rng_type;
extern gsl_rng* g_rng;
extern FILE* tvd_vs_gen_fstream; 

#define M_PI (3.14159265358979323846)
#define ONEOVERSQRT2PI  (1.0/sqrt(2.0*M_PI))

#define OUTPUT_TVD_VS_N (1)
#define OUTPUT_SAMPLES (1)
#define UNDERFLOW_ID (-1)
#define OVERFLOW_ID (-2)
#define TVD_EVERY (250)


// function declarations:
int propose_1dim(int i, int Width, int Ngrid_max);
double* propose(int Ndim, double* x_array, Proposal* prop, int is_ball);

double drand(void);

void print_array_of_int(int Nsize, int* array);
void print_array_of_double(int Nsize, double* array);

Ndim_histogram* init_target_distribution(int Ndim, int Ngrid_max, int normalize, const Binning_spec* bins);
double f_1dim (const Target_1dim* targ_1d, double x);
double F(const Target_1dim* targ_1d, int Ndim, double* x_array);
double draw_1dim(const Target_1dim* targ_1d);
double* draw_ndim(int n_dim, const Target_1dim* targ_1d);

