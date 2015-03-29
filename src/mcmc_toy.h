extern int n_accept;
extern int n_reject;
extern double sum_dsq;
extern Ndim_histogram* mcmc_out_hist;
extern int Ngrid_max;
extern double sigma;
extern int is_ball;
extern int n_modes;
extern double Xmin;
extern double Xmax;
extern Ndim_histogram* targp;
extern Targ_peak_1dim* peaks;
extern const gsl_rng_type* rng_type;
extern gsl_rng* the_rng;
FILE* tvd_vs_gen_fstream; 

#define OUTPUT_TVD_VS_N (1)
#define OUTPUT_SAMPLES (0)
#define UNDERFLOW_ID (-1)
#define OVERFLOW_ID (-2)

// function declarations:
double f_1dim (double x);
// double Fofint(int Ndim, int* index_array, int Ngrid_max);
double F(int Ndim, double* x_array);
int propose_1dim(int i, int Width, int Ngrid_max);
int* ipropose(int Ndim, int* index_array, Proposal* prop, int is_ball);
double* propose(int Ndim, double* x_array, Proposal* prop, int is_ball);
double drand(void);
void print_array_of_int(int Nsize, int* array);
void print_array_of_double(int Nsize, double* array);
int* i_array_from_x_array(int Ndim, double* x_array);
int i_from_x(double x);
double* x_array_from_i_array(int Ndim, int* i_array);
double x_from_i(int i);

int x_to_bin(Binning_spec* bin_spec, double x);

double* draw_ndim(int Ndim, int n_modes, Targ_peak_1dim* peaks);

double draw_1dim(int n_modes, Targ_peak_1dim* peaks);
