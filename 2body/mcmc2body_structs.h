// typedefs :
typedef struct{
  int n_dims;
  double* position;
  double width;
  double height;
  double shape; // shape param (bigger -> lighter tails)
}peak;

typedef struct{
  int n_peaks;
  peak** peaks;
}peaks;

typedef struct{
  int n_dims;
  double* xs;
}point;

typedef struct{
  int n_dims;
  int n_Ts;
  int n_walkers;
  double** w_xs;
  double* w_pis;
  int* w_ts; // 
  int* w_LRs;
  int* t_Lws;
  int* t_Rws;
  double* t_PIs;
  double* t_dsqrs;

  int* w_near_peak; // index of nearest peak
  int* w_transition_counts; // counts each walkers moves from one peak to the other.
  int* t_Lnear_peak;
  int* t_Ltransition_counts;
  int* t_Rnear_peak;
  int* t_Rtransition_counts;
  int* w_accepts; // counts the accepted X-changing moves for each walker
  int* t_accepts; // count the accepted X-changing moves for each T-level
  int* t_Laccepts;
  int* t_Raccepts;
  int* t_Tswap_accepts; // counts accepted T-swaps between T-levels i, i+1
int* t_Tswap_Laccepts;
int* t_Tswap_Raccepts;
}chain_state;

// function declarations:
//  peaks  function declarations
peaks* set_up_peaks(int n_dims, int n_peaks, double spacing, double width, double height_ratio, double shape_param);

//  chain_state  function declarations
chain_state* set_up_chain_state(int n_dims, int n_Ts, peaks* the_peaks, double init_width, double* inverse_Temperatures);
int check_state_consistency(chain_state* state, double* inverse_Temperatures);
void print_states_walker_order(chain_state* state); 
void print_states_T_order(chain_state* state);
void print_states_cold_only(chain_state* state);
void print_states_L0_Rtop_only(chain_state* state); // print L level 0, R level n_Ts-1

double pi(peaks* the_peaks, double* x);

void update_x(peaks* the_peaks, chain_state* state, double Tinverse, double prop_w, int iw);
void T_swap(chain_state* state, double* inverse_Temperatures);

void update_x_2b(peaks* the_peaks, chain_state* state, double Tinverse, double prop_w, int it);
void T_swap_2b_A(chain_state* state, double* inverse_Temperatures);
void T_swap_2b_B(chain_state* state, double* inverse_Temperatures);

double Dsquared(int n_dims, double* x, double* y);
double Kernel(double Dsq, double Lsq);
double PI(double pix, double piy, int it, int n_Ts, double K, double Tinverse);

// the end
