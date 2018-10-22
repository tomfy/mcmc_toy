mcmc2body: mcmc2body.o mcmc2body.h mcmc2body_structs.o mcmc2body_structs.h
	gcc -std=c99 -o mcmc2body mcmc2body.o mcmc2body_structs.o -lm -lgsl -lgslcblas

mcmc2body_structs.o: mcmc2body_structs.c mcmc2body_structs.h
	gcc -std=c99 -c mcmc2body_structs.c -O3

mcmc2body.o: mcmc2body.c mcmc2body.h mcmc2body_structs.h
	gcc -std=c99 -c mcmc2body.c -O3 


