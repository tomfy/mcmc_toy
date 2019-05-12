typedef struct{
  int n_dims;
  int type; // 0: gaussian, 1: 'cauchy'
  double* position;
  double* width; // 
  double height;
  double shape; // shape param (bigger -> lighter tails)
}peak;


typedef struct{
  int n_dims;
  int n_peaks;
  peak** peaks;
  //  long pi_evaluation_count;
  double* mean_x;
}target;


//  target  function declarations
target* set_up_target_from_argv(char** argv, int* i);
target* set_up_target_from_peaks(int n_peaks, peak** the_peaks);
target* set_up_target(int n_dims, int n_peaks, double spacing, double width, double height_ratio, double shape_param);

void print_target_info(const target* const the_target); // 
void print_peak_info(const peak* const a_peak); //
double min_peak_width(target* the_target);

double pi(const target* const the_target, const double* const x, int* which_peak);
