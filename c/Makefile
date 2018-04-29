#C makefile
all: diffusion newdiffusion

diffusion: diffusion.o
	$(CC) diffusion.o -o diffusion
  
diffusion.o: diffusion.c
	$(CC) $(CFLAGS) -std=c99 -O3 diffusion.c -c

newdiffusion: newdiffusion.c
	$(CC) $(CFLAGS) -std=c99 -O3 newdiffusion.c -lm -o newdiffusion
  
clean:
	rm *.o
	rm diffusion
	rm newdiffusion
