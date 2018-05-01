	/* CSC 435
	 * Final - Diffusion openmp
	 *
	 * Author: GnuForce
	 *
	 *
	 */

#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<omp.h>

	//declaration of the variables that define the system at the beginning
	//	and will be used throughout the program
	double mTotal, currTime, urms, D, rSize, rDiv, tStep, hval, conMax, conMin, hsqrd;

	//An int that will be used in a fashion similar to a boolean to
	//	control when the partition in activated
	int partition;

	//declaration of the step function that will perform a single iteration
	//	of the simulation
	void step (double* room, double* roomCopy, int* mask, int N);

	int main(){

		int nTemp = 0;
		char pValue;
		printf("How many divisions do want for the room?\n");
		scanf("%d",&nTemp);

		printf("Do you want a partition to be used?(y/n)\n");
		scanf(" %c",&pValue);

		int N = nTemp;
		double lRoom = 5;

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
		urms = 250.0;
		hval = (double) lRoom/N;
		D = 0.175;
		conMax = mTotal;
		conMin = 0.0;
		tStep = (double) (lRoom/urms)/N;
		hsqrd = (hval*hval);

		//declaration of for loop counter variables
		int i,j,k;

		//This represents the total molecules left after running the simulation
		//	used to check if matter consistency is held
		double tot = 0.0;

		//A 3 dimensional array that will operate as a rank 3 tensor used 
		//	to represent the room 
		double *room = calloc( 1, (N+2)*(N+2)*(N+2)*sizeof(double));
		double *roomCopy = calloc( 1, (N+2)*(N+2)*(N+2)*sizeof(double));
		int *mask = calloc( 1, (N+2)*(N+2)*(N+2)*sizeof(int));

		//Following for loops will initialize the room tensor with 0 values 
		//	when partioning is turned off, otherwise locations that 
		//	represent the partion in the room will be initialized
		//	to the value -1
		for (i=1; i<N+1; i++) {
			for (j=1; j<N+1; j++){
				for (k=1; k<N+1; k++){
		
					if(k == ((N+2)/2) && i >= ((N+2)/2)-1 && partition){
						mask[i*(2+N)*(2+N)+j*(2+N)+k] = 0;
					}else{
						mask[i*(2+N)*(2+N)+j*(2+N)+k] = 1;
					}

				}
			}
		}
	
		for( i=0;i<N+2;i++){
//			printf("\n");
			for(j =0;j<N+2;j++){
//				printf("\n");
				for(k=0;k<N+2;k++){
//					printf("%d ", mask[i*(2+N)*(2+N)+j*(2+N)+k]);
				}
			}
		}


		//Provides the room with the gas material to be dispersed
		//	to be understood as the "upper corner" of the room
		room[1*(2+N)*(2+N)+1*(2+N)+1] = mTotal;


		//We want the simulation to stop when the room has become sufficiently
		//	diffuse with the gas, thus we check if the ratio of lowest
		//	concentration to highest is less than 0.99, and when it is 
		//	higher we know the gas has diffused
		double* temp;
		while((conMin/conMax) < 0.99){
			currTime = currTime + tStep;

//			step(room, roomCopy, mask, N);



	//Every time we check to see the flux of gas between cells
	//	we would also need to multiply several values,
	//	slowing the speed of computation. By calculating the
	//	value once we only need to perform a single
	//	multiplication each time afterwards for each cell
	//	instead of several
	double coefficient = ((tStep*D) / (hsqrd) );


		omp_set_num_threads(1);

#pragma omp parallel shared(room,roomCopy,N) private(i,j,k)
{
	#pragma omp for
	for (i=1; i<N+1; i++){
		for (j=1; j<N+1; j++){

			for (k=1; k<N+1; k++){
				
				roomCopy[i*(2+N)*(2+N)+j*(2+N)+k] = room[i*(2+N)*(2+N)+j*(2+N)+k] + 
				(room[i*(2+N)*(2+N)+j*(2+N)+k+1] * mask[i*(2+N)*(2+N)+j*(2+N)+k+1] + room[i*(2+N)*(2+N)+j*(2+N)+k-1] * mask[i*(2+N)*(2+N)+j*(2+N)+k-1] + 
				room[i*(2+N)*(2+N)+j*(2+N)+k+(2+N)] * mask[i*(2+N)*(2+N)+j*(2+N)+k+(2+N)] + room[i*(2+N)*(2+N)+j*(2+N)+k-(2+N)] * mask[i*(2+N)*(2+N)+j*(2+N)+k-(2+N)] + 
				room[i*(2+N)*(2+N)+j*(2+N)+k+((2+N)*(2+N))] * mask[i*(2+N)*(2+N)+j*(2+N)+k+((2+N)*(2+N))] + room[i*(2+N)*(2+N)+j*(2+N)+k-((2+N)*(2+N))] * mask[i*(2+N)*(2+N)+j*(2+N)+k-((2+N)*(2+N))] -
				((mask[i*(2+N)*(2+N)+j*(2+N)+k+1] + mask[i*(2+N)*(2+N)+j*(2+N)+k-1] + mask[i*(2+N)*(2+N)+j*(2+N)+k+(2+N)] + mask[i*(2+N)*(2+N)+j*(2+N)+k-(2+N)] + mask[i*(2+N)*(2+N)+j*(2+N)+k+((2+N)*(2+N))] +  mask[i*(2+N)*(2+N)+j*(2+N)+k-((2+N)*(2+N))])
				 * room[i*(2+N)*(2+N)+j*(2+N)+k])) * 2 * coefficient;


			}
		}
	}
}


	//after resetting the concentration values we then find the values of min and max
	//	in order to tell when the loop shall end
	conMin = roomCopy[1*(2+N)*(2+N)+1*(2+N)+1];
	conMax = roomCopy[1*(2+N)*(2+N)+1*(2+N)+1];

	for (i=1; i<N+1; i++) {
		for (j=1; j<N+1; j++){
			for (k=1; k<N+1; k++){
				if (roomCopy[i*(2+N)*(2+N)+j*(2+N)+k] < conMin) {
					conMin = roomCopy[i*(2+N)*(2+N)+j*(2+N)+k];
				}
				if (roomCopy[i*(2+N)*(2+N)+j*(2+N)+k] > conMax) {
					conMax = roomCopy[i*(2+N)*(2+N)+j*(2+N)+k];
				}	
			}
		}
	}


			temp = room;
			room = roomCopy;
			roomCopy = temp;
		}

		//Here we total the values stored in all of the cells to check
		//	for any signifcant amount of lost or gained matter
		for (i=0; i<N+2; i++) {
			for (j=0; j<N+2; j++){
				for (k=0; k<N+2; k++){
					tot = tot + room[i*(2+N)*(2+N)+j*(2+N)+k];
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
		free(roomCopy);
		free(mask);
	}

