//extern int g_Ngrid_max;
//extern int g_is_ball;
extern const Target_1dim* g_targ_1d;
//extern const Ndim_histogram* g_targp;
extern const Ndim_histogram* g_targprobs_ndim;
extern const Ndim_histogram* g_targprobs_1dim;
extern const Ndim_histogram* g_targprobs_orthants;
extern const Ndim_histogram* g_targprobs_one_orthant;

extern const gsl_rng_type* g_rng_type;
extern gsl_rng* g_rng;
extern FILE* g_tvd_vs_gen_fstream; 

#define M_PI (3.14159265358979323846)
#define ONEOVERSQRT2PI  (1.0/sqrt(2.0*M_PI))

#define OUTPUT_SAMPLES (1)
#define UNDERFLOW_ID (-1)
#define OVERFLOW_ID (-2)
// #define OUTPUT_TVD_VS_N (1)
#define FIRST_SUMMARY_GEN (1000)
// #define TVD_EVERY (5000)
#define SUMMARY_GEN_FACTOR (1.3)
#define DO_EXACT (1)
#define DO_TVD (1)

// function declarations:
int propose_1dim(int i, int Width, int Ngrid_max);
double* propose(int Ndim, double* x_array, Proposal* prop);

double drand(void);

void print_array_of_int(int Nsize, int* array);
void print_array_of_double(int Nsize, double* array);

Ndim_histogram* init_target_distribution(int Ndim, int Ngrid_max, int normalize, const Binning_spec* bins);
double f_1dim (const Target_1dim* targ_1d, double x);
double F(const Target_1dim* targ_1d, int Ndim, double* x_array);
double draw_1dim(const Target_1dim* targ_1d);
double* draw_ndim(int n_dim, const Target_1dim* targ_1d);
double integral_f_1dim(const Target_1dim* targ_1d, double x, double y); // integral from x to y
double find_bin_upper_edge(const Target_1dim* targ_1d, double xlo, double Q);

double total_variation_distance(const Ndim_histogram* targprobs, const Ndim_histogram* hist);

double* merge_sorted_arrays(const int size1, const double* a1, const int size2, const double* a2);
double Kolmogorov_smirnov_D_statistic_2_sample(const int size1, const double* a1, const int size2, const double* a2);
double Kolmogorov_smirnov_D_statistic_1_sample(const int size1, const double* a1, double (*cdf)(double) );
double cdf(double y);
int cmpfunc (const void * a, const void * b);
double g(State* s);
double* copy_array(const int size, const double* a);

void add_arrays(int size, double* sum_x, const double* x);
