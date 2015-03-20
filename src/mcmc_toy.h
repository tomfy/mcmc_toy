// typedefs:
typedef int (*int_func_ptr)(); // don't need argument lists here?
typedef double (*double_func_ptr)();

typedef struct
{
  int n_dimensions;
  int* point; // point[0] = x index [0,Ngrid_max], point[1] = y index, etc.
  double prob; // 
  //  int xhalf; // 0 means  
  //  int yhalf; // 
} State;

typedef struct
{
  int W1; // 
  int W2;
  double p1; // prob of the W1 
} Proposal;

typedef struct
{
  State current_state;
  double temperature;
  int generation;
  Proposal proposal; // each T can have its own proposal
} Single_T_chain;

typedef struct
{
  int n_temperatures;
  Single_T_chain* coupled_chains; // array of metropolis-coupled chains at different T's
} Multi_T_chain;

typedef struct
{
  int n_replicates;
  Multi_T_chain* replicate_chains; 
} Replicate_chains;

typedef struct
{ // a Nsize^Ndim array
  int Ndim;
  int Nsize; 
  void* array; // an N dim array of N-1 dim arrays (or of double if N=1) 
} Ndim_array_of_double;

typedef struct
{ // a Nsize^Ndim array
  int Ndim;
  int Nsize; 
  void* array; // an N dim array of N-1 dim arrays (or of double if N=1) 
} Ndim_array_of_int;


typedef struct
{
  int Ndim; 
  int Ngrid_max; // grid [0,Ngrid_max] each dimension
  double total_weight;
  Ndim_array_of_double* weights; //
} Histogram_ndim;

/* typedef struct */
/* { */
/*   int Ndim;  */
/*   int Ngrid_max; // grid [0,Ngrid_max] each dimension */
/*   int sum; */
/*   Ndim_array_of_double* probs; // */
/* } Histogram_ndim; */

// function declarations:
// State
State* initialize_state(int n_dimensions, int Ngrid_max); //, Target_distribution* targp);
void free_state(State* s);
//int mcmc_step_2d(State* the_state, Proposal* prop, int Ngrid_max, Histogram_ndim* hist_ndim);

// Target_distribution* initialize_pi(int n_dimensions, int Ngrid_max);
//Histogram_2d* initialize_histogram_2d(int Ngrid_max);

// Histogram_ndim
Histogram_ndim* initialize_histogram_ndim(int Ndim, int Ngrid_max);
void normalize_pdf_ndim(Histogram_ndim* pdf);
double total_variation_distance(Histogram_ndim* targp, Histogram_ndim* mcmc_out);
Histogram_ndim* init_target_distribution(int Ndim, int Ngrid_max, int normalize);
void free_histogram_ndim(Histogram_ndim* h);
// double total_variation_distance(Target_distribution* targp, Histogram_ndim* hist_ndim);

// Ndim_array_of_double
Ndim_array_of_double* initialize_ndim_array_of_double(int Ndim, int Nsize, double init_value);
double* get_bin_pointer(Ndim_array_of_double* A, int* index_array);
double sum_ndim_array_of_double(Ndim_array_of_double* array_struct);
double set_ndim_array_of_double_with_function(Ndim_array_of_double* array_struct,
						     double outer_value, double_func_ptr function);
void multiply_ndim_array_of_double_by_scalar(Ndim_array_of_double* array_struct, double multiplier);
double sum_abs_difference_ndim_arrays_of_double(Ndim_array_of_double* a1, Ndim_array_of_double* a2);
void print_ndim_array_of_double(Ndim_array_of_double* A);
void free_ndim_array_of_double(Ndim_array_of_double* A);

// Ndim_array_of_int
Ndim_array_of_int* initialize_ndim_array_of_int(int Ndim, int Nsize, int init_value);


int mcmc_step_ndim(State* the_state, Proposal* prop, int Ngrid_max, Histogram_ndim* hist_ndim);
double f_1dim(double x);
double F(int Ndim, int* index_array, int Ngrid_max);
int propose_1dim(int i, int Width, int Ngrid_max);
int* propose(int Ndim, int* index_array, Proposal* prop);
double drand(void);
void print_int_array(int Nsize, int* array);

/* typedef struct */
/* { */
/*   double sum;  */
/*   int n_dimensions; // must be 2 for now. */
/*   int Ngrid_max; // grid [0,Ngrid_max] each dimension */
/*   Ndim_array_of_double* pi; // pi[i][j] is probability (unnormalized) of */
/* } Target_distribution; */
/* typedef struct */
/* { */
/*   int Ngrid_max; // grid [0,Ngrid_max] each dimension */
/*   int total_count; */
/*   int** bin_count; // */
/* } Histogram_2d; */
