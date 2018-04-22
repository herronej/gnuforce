//Code from Devin Mcbryde's 330 diffusion
//
//This code will use MPI on diffusion
#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>

#include<omp.h>
#include <mpi.h>

//#define N 100

//declaration of the variables that define the system at the beginning
//	and will be used throughout the program
double mTotal, currTime, mSpeed, D, rSize, rDiv, tStep, hval, conMax, conMin;

//An int that will be used in a fashion similar to a boolean to
//	control when the partition in activated
int partition;

// Put rank and size in global space
int rank, size;
int master;

//declaration of the step function that will perform a single iteration
//	of the simulation
void step (double* room, int N){
    	extern int rank, size; 
    	extern int MASTER;
    	MPI_Status *status;

    	// get the size of the parallel system and my rank
    	MPI_Comm_size(MPI_COMM_WORLD, &size);
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
}

int main()
{
  	int nTemp = 0;
	char pValue;
	printf("How many divisions do want for the room?\n");
	scanf("%d",&nTemp);

	printf("Do you want a partition to be used?(y/n)\n");
	scanf(" %c",&pValue);

	int N = nTemp;

	//printf("%c",pValue);

	if(pValue == 'y'){
		partition = 1;
	}else{
		partition = 0;
	}

	//cal for time objects that will tell us wall time at end of program
	time_t start;
	time_t end;
	//Give start the time of program starting execution
	start = time(NULL);

	//Initializing all of the necessary variables for the simulation to start
	mTotal = 1000000000000000000000.0;
	mSpeed = 250.0;
	hval = 5.0/N;
	D = 0.175;
	conMax = mTotal;
	conMin = 1.0;
	tStep = hval/mSpeed;
	
	//declaration of for loop counter variables
	int i,j,k;

	//This represents the total molecules left after running the simulation
	//	used to check if matter consistency is held
	double tot = 0.0;

	//A 3 dimensional array that will operate as a rank 3 tensor used 
	//	to represent the room 
	double *room = malloc(N*N*N*sizeof(double));

	//Following for loops will initialize the room tensor with 0 values 
	//	when partioning is turned off, otherwise locations that 
	//	represent the partion in the room will be initialized
	//	to the value -1
	for (i=0; i<N; i++) {
		for (j=0; j<N; j++){
			for (k=0; k<N; k++){
				//negative values are used when partition is true and will
				//	place them half way into the room (when j == (N/2)-1)
				//	and half way up (when i >= (N/2)-1)
				if(j == (N/2)-1 && i >= (N/2)-1 && partition){
					room[i*N*N+j*N+k] = -1.0;
				}else{
					room[i*N*N+j*N+k] = 0.0;
				}
			}
		}
	}

  	//Provides the room with the gas material to be dispersed
	//	to be understood as the "upper corner" of the room
	room[0] = mTotal;

	//We want the simulation to stop when the room has become sufficiently
	//	diffuse with the gas, thus we check if the ratio of lowest
	//	concentration to highest is less than 0.99, and when it is 
	//	higher we know the gas has diffused
	while((conMin/conMax) < 0.99){
		currTime = currTime + tStep;
		step(room, N);
	}

	//Here we total the values stored in all of the cells to check
	//	for any signifcant amount of lost or gained matter
	for (i=0; i<N; i++) {
                for (j=0; j<N; j++){
                        for (k=0; k<N; k++){
                                tot = tot + room[i*N*N+j*N+k];
                        }
                }
        }


	//output of the simulation detailing 5 vaules
	//	How many molecules did we start with
	//	How many molecules did we end with
	//	The total amount of time it took for the room to become diffused
	//	The minimum concentration in the room
	//	the maximum concentration in the room
	
	end = time(NULL);
	double seconds = difftime(end,start);
	printf("Total molecules starting: %f\n", mTotal);
	printf("Total molecules left: %f\n", tot);
	printf("Time Simulated: %f\n", currTime);
	printf("min concentration: %f\n", conMin);
	printf("max concentration: %f\n", conMax);
	printf("Wall time: %f\n", seconds);
	//free the memory held by the array
	free(room);

	return 0;
}
