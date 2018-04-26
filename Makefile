#Top level makefile
all : fortran c

#fdiff: fdiff.o 
#	$(F90) fdiff.o -o fdiff

#fdiff.o: fdiff.f90
#	$(F90) fdiff.f90 -c

#cdiff: cdiff.o
#	$(CC) cdiff.o -o cdiff

#cdff.o: cdiff.c
#	$(CC) cdiff.c -c 

c:
	cd c && $(MAKE)

fortran:
	cd fortran && $(MAKE)

clean :
	cd fortran && $(MAKE) clean
	cd c && $(MAKE) clean
	#rm *.o

#This next target get "made" every time
.PHONY: fortran c
