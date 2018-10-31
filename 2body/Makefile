mcmc2body: mcmc2body.o mcmc2body.h mcmc2body_structs.o mcmc2body_structs.h target.o target.h chain_architecture.o chain_architecture.h
	gcc -std=c99 -o mcmc2body mcmc2body.o target.o chain_architecture.o mcmc2body_structs.o -lm -lgsl -lgslcblas

mcmc2body_structs.o: mcmc2body_structs.c mcmc2body_structs.h target.h chain_architecture.h
	gcc -std=c99 -c mcmc2body_structs.c -O3

chain_architecture.o: chain_architecture.c chain_architecture.h target.h mcmc2body_structs.h
	gcc -std=c99 -c chain_architecture.c -O3

target.o: target.h target.h 
	gcc -std=c99 -c target.c -O3

mcmc2body.o: mcmc2body.c mcmc2body.h target.h mcmc2body_structs.h chain_architecture.h
	gcc -std=c99 -c mcmc2body.c -O3 


