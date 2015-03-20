// function declarations:

int mcmc_step_ndim(State* the_state, Proposal* prop, int Ngrid_max, Ndim_histogram* hist_ndim);
double f_1dim(double x);
double F(int Ndim, int* index_array, int Ngrid_max);
int propose_1dim(int i, int Width, int Ngrid_max);
int* propose(int Ndim, int* index_array, Proposal* prop, int is_ball);
double drand(void);
void print_int_array(int Nsize, int* array);
