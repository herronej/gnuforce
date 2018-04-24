all : cdiff fdiff

fdiff: fdiff.o 
	$(F90) mpif90 -std=f90 fdiff.o -o fdiff

fdiff.o: fdiff.f90
	$(F90) fdiff.f90 -c

cdiff: cdiff.o
	$(CC) mpicc -std=c99 cdiff.o -o cdiff

#cdiff: cdiff.o
#        $(CC) mpicc -std=c99 cdiff.o -o cdiff -lm

cdff.o: cdiff.c
	$(CC) cdiff.c -c 

clean :
	rm *.o
	rm *.a
