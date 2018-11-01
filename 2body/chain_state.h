// typedefs :

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
chain_state* set_up_chain_state(int n_dims, int n_Ts, target* the_target, chain_architecture* arch, double init_width);
int check_state_consistency(target* the_target, chain_architecture* arch, chain_state* state);
void print_states_walker_order(chain_state* state); 
void print_states_T_order(chain_state* state);
void print_states_cold_only(chain_state* state);
void print_states_L0_Rtop_only(chain_state* state); // print L level 0, R level n_Ts-1

// ***** updates
// 1-body:
void update_x(target* the_target, chain_architecture* arch, chain_state* state, int iw);
void T_swap(chain_state* state, double* inverse_Temperatures);


// 2-body
void LR_swap(chain_architecture* arch, chain_state* state);
void update_x_2b(target* the_target, chain_architecture* arch, chain_state* state, int it);
void T_swap_2b_A(chain_architecture* arch, chain_state* state);
void T_swap_2b_B(chain_architecture* arch, chain_state* state);

void step_1b(target* target, chain_architecture* arch, chain_state* state);
void step_2b(target* target, chain_architecture* arch, chain_state* state);

void cold_transition_observe_and_count_sym(chain_state* state);
void cold_transition_observe_and_count_asym(chain_state* state);



// the end
