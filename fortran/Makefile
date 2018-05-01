#Fortran makefile
F90=gfortran
FFLAGS= -O3 -fopenmp
all : cube_mem.o omp_diffusion
#serial: serial_diffusion.o 
#	$(F90) serial_diffusion.o -o serial 
omp_diffusion : omp_diffusion.f90
	$(F90) $(FFLAGS) cube_mem.o omp_diffusion.f90 -o omp_diffusion

cube_mem.o : cube_mem.f90
	$(F90) $(FFLAGS) -c cube_mem.f90
 
clean:
	rm *.o
	rm *.mod
	rm omp_diffusion
