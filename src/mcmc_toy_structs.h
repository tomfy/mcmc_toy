// typedefs:
typedef int (*int_func_ptr)(); // don't need argument lists here?
typedef double (*double_func_ptr)();


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
} Ndim_histogram;

typedef struct {
  double position;
  double sigma;
  double weight;
} Targ_peak_1dim;

// mcmc specific structs
typedef struct
{
  int n_dimensions;
  int* ipoint; // ipoint[0] = x index [0,Ngrid_max], ipoint[1] = y index, etc.
  double* point; //
  double prob; // 
} State;

typedef struct
{
  double W1; // 
  double W2;
  double p1; // prob of the W1 
} Proposal;

typedef struct
{
  State* current_state;
  double temperature;
  int generation;
  Proposal proposal; // each T can have its own proposal
  Ndim_histogram* mcmc_out_hist;
} Single_T_chain;

typedef struct
{
  int n_temperatures;
  Single_T_chain** coupled_chains; // array of metropolis-coupled chains at different T's
} Multi_T_chain;

typedef struct
{
  int n_replicates;
  Multi_T_chain* replicate_chains; 
} Replicate_chains;

typedef struct
{    // < x_lo -> 'underflow'; > x_hi -> 'overflow' x_lo <= x < bin_ubs[0] -> bin 0, 
  int n_bins; 
  double x_lo;
  double* bin_edges; // bin_edges[0] is lower edge of bin 0, etc. bin_edges[i] is upper edge of bin i-1
  double x_hi;
} Binning_spec;

// function declarations:

// State
State* construct_state(int n_dimensions, int Ngrid_max); //, Target_distribution* targp);
void free_state(State* s);

// Single_T_chain
Single_T_chain* construct_single_T_chain(double T, Proposal P, State* S);
State* single_T_chain_mcmc_step(Single_T_chain* chain);
void free_single_T_chain(Single_T_chain* chain);
// Multi_T_chain
Multi_T_chain* construct_multi_T_chain(int n_temperatures, double* temperatures, Proposal* proposals, State** states);
void multi_T_chain_within_T_mcmc_step(Multi_T_chain* multi_T_chain);
void multi_T_chain_T_swap_mcmc_step(Multi_T_chain* multi_T_chain);
void free_multi_T_chain(Multi_T_chain* chain);
// Ndim_histogram
Ndim_histogram* construct_ndim_histogram(int Ndim, int Ngrid_max);
Ndim_histogram* construct_copy_ndim_histogram(Ndim_histogram* A);
void normalize_ndim_histogram(Ndim_histogram* pdf);
double total_variation_distance(Ndim_histogram* targp, Ndim_histogram* mcmc_out);
Ndim_histogram* init_target_distribution(int Ndim, int Ngrid_max, int normalize);
void* free_ndim_histogram(Ndim_histogram* h);

// Ndim_array_of_double
Ndim_array_of_double* construct_ndim_array_of_double(int Ndim, int Nsize, double init_value);
Ndim_array_of_double* construct_copy_ndim_array_of_double(Ndim_array_of_double* A);
double* get_pointer_to_element(Ndim_array_of_double* A, int* index_array);
double sum_ndim_array_of_double(Ndim_array_of_double* array_struct);
double set_ndim_array_of_double_with_function(Ndim_array_of_double* array_struct,
						     double outer_value, double_func_ptr function);
double set_ndim_array_of_double_to_target(Ndim_array_of_double* array_struct,
					  double outer_value, Targ_peak_1dim* peaks, Binning_spec* bins);
void multiply_ndim_array_of_double_by_scalar(Ndim_array_of_double* array_struct, double multiplier);
double sum_abs_difference_ndim_arrays_of_double(Ndim_array_of_double* a1, Ndim_array_of_double* a2);
void print_ndim_array_of_double(Ndim_array_of_double* A);
void free_ndim_array_of_double(Ndim_array_of_double* A);

// Ndim_array_of_int
Ndim_array_of_int* construct_ndim_array_of_int(int Ndim, int Nsize, int init_value);

// array of Targ_peak_1dim 
void normalize_targ_1dim(int n_modes, Targ_peak_1dim* peaks);

//
Binning_spec* construct_binning_spec(int n_bins, Targ_peak_1dim* peaks);
