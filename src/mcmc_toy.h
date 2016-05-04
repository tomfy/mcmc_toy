extern const Target_1dim* g_targ_1d;
extern const Ndim_histogram* g_targprobs_ndim;
extern const Ndim_histogram* g_targprobs_1dim;
extern const Ndim_histogram* g_targprobs_orthants;
extern const Ndim_histogram* g_targprobs_one_orthant;

extern const gsl_rng_type* g_rng_type;
extern gsl_rng* g_rng;
extern FILE* g_tvd_vs_gen_fstream; 
extern FILE* g_run_params_fstream;
extern FILE* g_accept_info_fstream;

extern long g_n_pi_evaluations;
extern long g_max_pi_evaluations;


#define M_PI (3.14159265358979323846)
#define SQRT2PI  (sqrt(2.0*M_PI))
#define ONEOVERSQRT2PI  (1.0/sqrt(2.0*M_PI))

#define OUTPUT_SAMPLES (0)
#define FIRST_SUMMARY_GEN (1000)
#define SUMMARY_GEN_FACTOR (2.0)

#define USELOGPROB (1)
// #define DO_TVD (1)

// function declarations:
int propose_1dim(int i, int Width, int Ngrid_max);
// double* propose(int Ndim, double* x_array, Proposal* prop);
double* propose(int n_dim, double* x_array, Proposal* prop, int* which_prop);
double drand(void);

void print_array_of_int(int Nsize, int* array);
void print_array_of_double(int Nsize, double* array);

Ndim_histogram* init_target_distribution(int Ndim, int Ngrid_max, int normalize, const Binning_spec* bins);
double f_1dim (const Target_1dim* targ_1d, double x);
double f_ndim(const Target_1dim* targ_1d, int Ndim, double* x_array);
double log_f_ndim(const Target_1dim* targ_1d, int n_dim, double* x_array);
double log_f_1dim(const Target_1dim* targ_1d, double x);
double draw_1dim(const Target_1dim* targ_1d, double temperature);
double* draw_ndim(int n_dim, const Target_1dim* targ_1d, double temperature);
double integral_f_1dim(const Target_1dim* targ_1d, double x, double y); // integral from x to y
double cdf(double y);

double find_bin_upper_edge(const Target_1dim* targ_1d, double xlo, double Q);

double total_variation_distance(const Ndim_histogram* targprobs, const Ndim_histogram* hist);

double* merge_sorted_arrays(const int size1, const double* a1, const int size2, const double* a2);
double Kolmogorov_smirnov_D_statistic_2_sample(const int size1, const double* a1, const int size2, const double* a2);
double Kolmogorov_smirnov_D_statistic_1_sample(const int size1, const double* a1, double (*cdf)(double) );
double edf_statistic(const int size, const double* a, double (*cdf)(double) );
double anderson_darling_statistic1(const int size, const double* a, double (*cdf)(double) );
double anderson_darling_statistic(const int size, const double* a, double (*cdf)(double) );

int cmpfunc (const void * a, const void * b);
double g(const State* s);
double g_diag(const State* s);
double g_shortrange(const State* s);
double* copy_array(const int size, const double* a);

void add_arrays(int size, double* sum_x, const double* x);
