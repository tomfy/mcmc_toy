//

extern const gsl_rng_type* g_rng_type;
extern gsl_rng* g_rng;

extern double two_body_interpolation_power; // 0 -> geom, 1 -> linear, other -> weighted power mean with power P.
extern long pi_evaluation_count;

double Dsquared(int n_dims, double* x, double* y);

// the end
