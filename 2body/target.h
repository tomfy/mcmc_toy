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
}target;


//  target  function declarations
target* set_up_target(int n_dims, int n_peaks, double spacing, double width, double height_ratio, double shape_param);
void print_target_info(target* the_target); // not implemented.
double pi(target* the_target, double* x);
