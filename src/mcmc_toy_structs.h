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
{    // < x_lo -> 'underflow'; > x_hi -> 'overflow' x_lo <= x < bin_ubs[0] -> bin 0, 
  int n_bins; 
  double x_lo;
  double* bin_edges; // bin_edges[0] is lower edge of bin 0, etc. bin_edges[i] is upper edge of bin i-1
  double x_hi;
} Binning_spec;

typedef struct
{
  Binning_spec* ndim_bins;
  Binning_spec* onedim_bins; // for n_dim > 2, do n_dim 2d histograms
  Binning_spec* orthants_bins; // 1 bin for each orthant (i.e. n_bin = 2, each dimension if weights equal)
  Binning_spec* positive_bins; //
} Binning_spec_set;

typedef struct
{
  int Ndim;
  const Binning_spec* bins; // specifies where bin edges are (same for all dimensions)
  int Ngrid_max; // grid [0,Ngrid_max] each dimension
  double total_weight;
  Ndim_array_of_double* weights; //
} Ndim_histogram;

typedef struct {
  double position;
  double sigma;
  double weight;
} Target_peak_1dim;

typedef struct {
  int n_modes;
  Target_peak_1dim* peaks;
} Target_1dim;

// mcmc specific structs
typedef struct
{
  int n_dimensions;
  //  int* ipoint; // ipoint[0] = x index [0,Ngrid_max], ipoint[1] = y index, etc.
  double* point; //
  double prob; // 
} State;

typedef struct
{
  char* shape; // "cube", "ball", "gaussian"
  double W1; // 
  double W2;
  double p1; // prob of the W1 
} Proposal;

typedef struct
{
  State* current_state;
  double temperature;
  int generation;
  int n_accept;
  int n_reject;
  double dsq_sum;
  int n_jumps;
  Proposal proposal; // each T can have its own proposal
  Ndim_histogram* mcmc_out_hist;
  Ndim_histogram** mcmc_out_1d_hists; // for n_dim > 2, do n_dim 2d histograms
  Ndim_histogram* mcmc_out_orthants_hist; // 1 bin for each orthant (i.e. n_bin = 2, each dimension if weights equal)
  Ndim_histogram* mcmc_out_reflected_hist; // pts reflected into all-positive orthant

  Ndim_histogram* exact_draw_hist;
  Ndim_histogram** exact_draw_1d_hists; // for n_dim > 2, do n_dim 2d histograms
  Ndim_histogram* exact_draw_orthants_hist; // 1 bin for each orthant (i.e. n_bin = 2, each dimension if weights equal)
  Ndim_histogram* exact_draw_reflected_hist; // pts reflected into all-positive orthant
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

// function declarations:

// State
//State* construct_state(int n_dimensions); //, Target_distribution* targp);
State* construct_state(int n_dimensions, const Target_1dim* targ_1d); 
void free_state(State* s);

// Single_T_chain
Single_T_chain* construct_single_T_chain(double T, Proposal P, State* S, const Binning_spec_set* bins);
State* single_T_chain_mcmc_step(Single_T_chain* chain);
void single_T_chain_histogram_current_state(Single_T_chain* chain);
void single_T_chain_output_tvd(Single_T_chain* chain);
void free_single_T_chain(Single_T_chain* chain);
// Multi_T_chain
Multi_T_chain* construct_multi_T_chain(int n_temperatures, double* temperatures, Proposal* proposals, State** states, const Binning_spec_set* bins);
void multi_T_chain_within_T_mcmc_step(Multi_T_chain* multi_T_chain);
void multi_T_chain_T_swap_mcmc_step(Multi_T_chain* multi_T_chain);
void free_multi_T_chain(Multi_T_chain* chain);
// Ndim_histogram
Ndim_histogram* construct_ndim_histogram(int n_dim, const Binning_spec* bins); //(int Ndim, int Ngrid_max);
Ndim_histogram* construct_copy_ndim_histogram(const Ndim_histogram* A);
void normalize_ndim_histogram(Ndim_histogram* pdf);
// double total_variation_distance(const Ndim_histogram* targp, const Ndim_histogram* mcmc_out);
//Ndim_histogram* init_target_distribution(int Ndim, int Ngrid_max, int normalize);

void add_data_pt_to_ndim_histogram(Ndim_histogram* A, int n_dim, double* xs);
void* free_ndim_histogram(Ndim_histogram* h);

// Ndim_array_of_double
Ndim_array_of_double* construct_ndim_array_of_double(int Ndim, int Nsize, double init_value);
Ndim_array_of_double* construct_copy_ndim_array_of_double(Ndim_array_of_double* A);
double* get_pointer_to_element(Ndim_array_of_double* A, int* index_array);
double sum_ndim_array_of_double(Ndim_array_of_double* array_struct);
double set_ndim_array_of_double_with_function(Ndim_array_of_double* array_struct,
						     double outer_value, double_func_ptr function);
double set_ndim_array_of_double_to_target(Ndim_array_of_double* array_struct,
					  double outer_value, Target_peak_1dim* peaks, const Binning_spec* bins);
void multiply_ndim_array_of_double_by_scalar(Ndim_array_of_double* array_struct, double multiplier);
double sum_abs_difference_ndim_arrays_of_double(const Ndim_array_of_double* a1, const Ndim_array_of_double* a2);
void print_ndim_array_of_double(const Ndim_array_of_double* A);
void free_ndim_array_of_double(Ndim_array_of_double* A);

// Ndim_array_of_int
Ndim_array_of_int* construct_ndim_array_of_int(int Ndim, int Nsize, int init_value);

// Target_1dim 
void normalize_targ_1dim(Target_1dim* targ);

// Binning_spec
Binning_spec* construct_binning_spec_old(int n_bins, const Target_1dim* targ_1d);
Binning_spec* construct_binning_spec(int n_bins, const Target_1dim* targ_1d, double xlo, double xhi);
int x_to_bin(const Binning_spec* bin_spec, double x);
int* i_array_from_x_array(const Binning_spec* bins, int Ndim, double* x_array);


void print_tvd(Ndim_histogram* hist, const Ndim_histogram* targprobs);
