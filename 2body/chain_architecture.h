typedef struct{
  int n_per_level; // 1: 'classic' heating; 2: 2-body;
  int symmetry; // 1: symmetrical, 0: asymmetrical.
  int n_levels; // = n_Ts for standard **(1/T) heating.
  double two_body_interpolation_power;
  double* inverse_Temperatures;
  double* Lprop_widths;
  double* Rprop_widths;
  double* kernel_widths;
}chain_architecture;

chain_architecture* set_up_chain_architecture(int n_levels, int n_per_level, int symmetry, double Thot, double min_prop_w, double max_prop_w, double min_kernel_width, double max_kernel_width, double two_body_interpolation_power);
void print_chain_architecture_info(const chain_architecture* const arch);

double Kernel(double Dsq, double Lsq);
double PI(double pix, double piy, int it, int n_Ts, double K, const chain_architecture* const arch);
