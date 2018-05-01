#C makefile
CC=gcc
CFLAGS= -fopenmp -std=c99 -O3 
all: diffusion

diffusion: diffusion.c
	$(CC) $(CFLAGS) -o diffusion  diffusion.c

#diffusion: diffusion.o
#	$(CC) diffusion.o -o diffusion

  
clean:
	rm diffusion
	rm *.o
