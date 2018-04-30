#C makefile
all: diffusion

diffusion: diffusion.o
	$(CC) diffusion.o -o diffusion
  
diffusion.o: diffusion.c
	$(CC) $(CFLAGS) -std=c99 -O3 diffusion.c -c

#diffusion: diffusion.o
#	$(CC) diffusion.o -o diffusion

  
clean:
	rm *.o
	rm diffusion
