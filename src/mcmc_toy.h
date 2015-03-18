// typedefs:
typedef struct
{
  int n_dimensions;
  int* point; // point[0] = x index [0,Ngrid_max], point[1] = y index, etc.
  double prob; // 
  int xhalf; // 0 means  
  int yhalf; // 
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
{
  int n_dimensions; // must be 2 for now.
  int Ngrid_max; // grid [0,Ngrid_max] each dimension
  double** pi; // pi[i][j] is probability (unnormalized) of point i,j
} Target_distribution;

typedef struct
{
  int Ngrid_max; // grid [0,Ngrid_max] each dimension
  int total_count;
  int** bin_counts; //
} Histogram_2d;

// function declarations:
State* initialize_state(int n_dimensions, int Ngrid_max, Target_distribution* targp);
int mcmc_step(State* the_state, Proposal* prop, Target_distribution* targp, Histogram_2d* hist2d);
double f(double x);
int propose_1dim(int i, int Width, int Ngrid_max);
double drand(void);
Target_distribution* initialize_pi(int n_dimensions, int Ngrid_max);
Histogram_2d* initialize_histogram_2d(int Ngrid_max);
double total_variation_distance(Target_distribution* targp, Histogram_2d* hist2d);
