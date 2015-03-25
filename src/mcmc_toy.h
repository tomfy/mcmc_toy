extern int n_accept;
extern int n_reject;
extern double sum_dsq;
extern Ndim_histogram* mcmc_out_hist;
extern int Ngrid_max;
extern double sigma;
extern int is_ball;
extern int n_modes;


// function declarations:

//int mcmc_step(State* the_state, Proposal* prop, int Ngrid_max, Ndim_histogram* hist_ndim);
double f_1dim(double x);
double F(int Ndim, int* index_array, int Ngrid_max);
int propose_1dim(int i, int Width, int Ngrid_max);
int* propose(int Ndim, int* index_array, Proposal* prop, int is_ball);
double drand(void);
void print_array_of_int(int Nsize, int* array);
