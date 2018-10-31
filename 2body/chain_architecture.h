typedef struct{
  int n_per_level; // 1: 'classic' heating; 2: 2-body;
  int symmetry; // 1: symmetrical, 0: asymmetrical.
  int n_levels; // = n_Ts for standard **(1/T) heating.
  double* inverse_Temperatures;
  double* Lprop_widths;
  double* Rprop_widths;
  double* kernel_widths;
}chain_architecture;

chain_architecture* set_up_chain_architecture(int n_per_level, int symmetry, int n_levels, double Thot, double min_prop_w, double max_prop_w, double min_kernel_width, double max_kernel_width);
void print_chain_architecture_info(chain_architecture* arch);

