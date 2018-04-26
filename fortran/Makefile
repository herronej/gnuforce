#Fortran makefile
all : cube_mem.o serial_diffusion

#serial: serial_diffusion.o 
#	$(F90) serial_diffusion.o -o serial 

cube_mem.o : cube_mem.f90
	$(F90) $(FFLAGS) gfortran -c cube_mem.f90
 
serial_diffusion : serial_diffusion.f90
	$(F90) $(FFLAGS) gfortran cube_mem.o serial_diffusion.f90 -o serial_diffusion
 
clean:
	rm *.o
	rm *.mod
	rm serial_diffusion
