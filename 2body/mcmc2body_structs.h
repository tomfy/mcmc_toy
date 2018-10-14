// typedefs :
typedef struct{
  int n_dims;
  double* position;
  double width;
  double height;
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
}chain_state;

// function declarations:


void old_print_state_walker_order(int n_Ts, int n_d, double** xs); 
void print_states_T_order(chain_state* state);
double pi(peaks* the_peaks, double* x);
void update_x(int n_dims, peaks* the_peaks,  int iw, chain_state* state, double Tinverse);
void T_swap(chain_state* state, double* inverse_Temperatures);
int check_state_consistency(int n_dims, int n_Ts, double* pis, double** xs,  int* Ts, int* LRs, int* Lws, int* Rws);


// the end
