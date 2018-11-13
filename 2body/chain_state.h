// typedefs :

typedef struct{
  int n_dims;
  double* xs;
}point;

typedef struct{
  int n_dims;
  int n_levels;
  int n_walkers;
  double** w_xs;
  double* w_pis;
  int* w_ts; // 
  int* w_LRs;
  int* t_Lws;
  int* t_Rws;
  double* t_PIs;
  double* t_dsqrs;

  long updates; // the number of updates so far.
  double** t_Lxsums;
  double** t_Rxsums;
  vector* all_coldx0s;
  vector* all_coldx1s;
  vector** w_coldx0s;
  vector** w_coldx1s; 
  int* w_near_peak; // index of nearest peak
  int* w_transition_counts; // counts each walker's moves from one peak to the other.
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
  int* t_LRswap_accepts; // counts accepted LR-swaps for each T-level. 
}chain_state;

// function declarations:

//  chain_state  function declarations
chain_state* set_up_chain_state(const target* const the_target, const chain_architecture* const arch, double init_width);
void reset_chain_state(chain_state* state); // reset the cumulative variables to initial values
void free_chain_state(chain_state* state);

int check_state_consistency(const target* const the_target, const chain_architecture* const arch, const chain_state* const state);
void print_states_walker_order(const chain_state* const state); 
void print_states_T_order(const chain_state* const state);
void print_states_cold_only(const chain_state* const state);
void print_states_L0_Rtop_only(const chain_state* const state); // print L level 0, R level n_Ts-1

// ***** updates
// 1-body:
void update_x(const target* const the_target, const chain_architecture* const arch, chain_state* state, int iw);
void T_swap(chain_state* state, double* inverse_Temperatures);


// 2-body
void LR_swap(const chain_architecture* const arch, chain_state* state);
void update_x_2b(const target* const the_target, const chain_architecture* const arch, chain_state* state, int it);
void T_swap_2b_A(const chain_architecture* const arch, chain_state* state);
void T_swap_2b_B(const chain_architecture* const arch, chain_state* state);

void step_1b(const target* const target, const chain_architecture* const arch, chain_state* state);
void step_2b(const target* const target, const chain_architecture* const arch, chain_state* state);

void cold_transition_observe_and_count_sym(chain_state* state);
void cold_transition_observe_and_count_asym(chain_state* state);

void accumulate_x_sums(chain_state* state);


// the end
