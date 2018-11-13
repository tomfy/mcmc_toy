/* typedef struct */
/* { // a Nsize^Ndim array */
/*   int Ndim; */
/*   int Nsize;  */
/*   void* array; // an N dim array of N-1 dim arrays (or of double if N=1)  */
/* } Ndim_array_of_double; */


/* typedef struct */
/* {    // < x_lo -> 'underflow'; > x_hi -> 'overflow' x_lo <= x < bin_ubs[0] -> bin 0,  */
/*   int n_bins;  */
/*   double x_lo; */
/*   double* bin_edges; // bin_edges[0] is lower edge of bin 0, etc. bin_edges[i] is upper edge of bin i-1 */
/*   double x_hi; */
/* } Binning_spec; */

/* typedef struct */
/* { */
/*   int Ndim; */
/*   const Binning_spec* bins; // specifies where bin edges are (same for all dimensions) */
/*   int Ngrid_max; // grid [0,Ngrid_max] each dimension */
/*   double total_weight; */
/*   Ndim_array_of_double* weights; // */
/* } Ndim_histogram; */

/* // Ndim_histogram */
/* Ndim_histogram* construct_ndim_histogram(int n_dim, const Binning_spec* bins); //(int Ndim, int Ngrid_max); */
/* Ndim_histogram* construct_copy_ndim_histogram(const Ndim_histogram* A); */
/* void normalize_ndim_histogram(Ndim_histogram* pdf); */


// ********************************************************************


typedef struct{
  long capacity;
  long count;
  double* elements;
}vector;

vector* new_vector(long capacity);
int push(vector* v, double x);
int pop(vector* v, double* x);
void free_vector(vector* v);


double Dsquared(int n_dims, double* x, double* y);
double error(int n_dims, double* mean_x, double* sum_x, long n);
double* sum_arrays(int n_dims, double* x1, double* x2); // get a new array which is the element by element sum of two arrays (of same size)
double* add_array_to_array(int n_dims, double* x1, double* x2); // to each element of x1, add the corresponding elem. of x2
double* mult_array_by_scalar(double A, int n_dims, double* x); // multiply each element of x by A

int compare_doubles(const void * a, const void * b);


double Kolmogorov_Smirnov_D_statistic_2_sample(const int size1, const double* a1, const int size2, const double* a2);
