#C makefile
CC=gcc
CFLAGS= -fopenmp -std=c99 -O3 -march=native -mtune=native
all: diffusion

diffusion: diffusion.c
	$(CC) $(CFLAGS) -o diffusion  diffusion.c

#diffusion: diffusion.o
#	$(CC) diffusion.o -o diffusion

  
clean:
	rm diffusion
	rm *.o
