# mcmc_toy
MCMC and MCMCMC. Space is [0..N]^d , proposals are symmetric. 

compile with gcc  -std=c99  mcmc_toy.c -lm

Target density is presently of form Prod(f(x_i)), i.e. there is a single 1-d function f, and
the d-dimensionl prob. is  f(x1)*f(x2)* ... * f(x_d).

Proposal is presently uniform over the ((2W+1)^d - 1) pts in d-dim cube centered on present pt, 
or mixture of two such proposals with different widths.

The idea is to test things like what is best proposal width (or mixture of widths) for various 
target distributions, particularly multimodal.

