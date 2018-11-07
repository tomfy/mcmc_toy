mcmc2body: mcmc2body.o mcmc2body.h target.o target.h chain_architecture.o chain_architecture.h chain_state.o chain_state.h globals.h util.o util.h
	gcc -std=c99 -o mcmc2body mcmc2body.o target.o chain_architecture.o chain_state.o util.o -lm -lgsl -lgslcblas

chain_state.o: chain_state.c chain_state.h target.h chain_architecture.h globals.h util.o util.h
	gcc -std=c99 -c chain_state.c -O3

chain_architecture.o: chain_architecture.c chain_architecture.h target.h
	gcc -std=c99 -c chain_architecture.c -O3

target.o: target.c target.h globals.h util.o util.h
	gcc -std=c99 -c target.c -O3

util.o: util.c util.h
	gcc -std=c99 -c util.c -O3

mcmc2body.o: mcmc2body.c mcmc2body.h target.h chain_architecture.h chain_state.h globals.h util.o util.h 
	gcc -std=c99 -c mcmc2body.c -O3 


