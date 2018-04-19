/* CSC 330
 * Assignment 2 - Diffusion
 *
 * Author: Devin McBryde
 *
 *
 */

#include<stdlib.h>
#include<stdio.h>
#include<time.h>

//declaration of the variables that define the system at the beginning
//	and will be used throughout the program
double mTotal, currTime, mSpeed, D, rSize, rDiv, tStep, hval, conMax, conMin;

//An int that will be used in a fashion similar to a boolean to
//	control when the partition in activated
int partition;

//declaration of the step function that will perform a single iteration
//	of the simulation
void step (double* room, int N);

int main(){

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
}

void step(double* room, int N){

	//The vaules used to iterate through the room array
	//	using nested for loops
	int i,j,k;
	
	//This array will store the different values of 
	//	change in concentration between two cells
	//	The name means concentration difference
	double dCon[6];
	for(i=0;i<6;i++){dCon[i] = 0.0;}

	//Every time we check to see the flux of gas between cells
	//	we would also need to multiply several values,
	//	slowing the speed of computation. By calculating the
	//	value once we only need to perform a single
	//	multiplication each time afterwards for each cell
	//	instead of several
	double coefficient = ((tStep*D) / (hval*hval));

	for (i=0; i<N; i++){
                for (j=0; j<N; j++){
			for (k=0; k<N; k++){


				//calculate the difference in concentration from flux with each cube face
				//	The 6 faces of the cube are represented with different address values
				//	and an if is used to determine if it is safe to move molecules
				//	if a value = N-1  or 0 then we have hit a face of the cube and do not calculate
				if(room[i*N*N+j*N+k] != -1){
				if(k==N-1 || room[i*N*N+j*N+k+1] == -1){
					dCon[0] = 0;
				}else{
					dCon[0] = (room[i*N*N+j*N+k]    -room[i*N*N+j*N+k+1]) * coefficient;
					room[i*N*N+j*N+k] = room[i*N*N+j*N+k] - dCon[0];
					room[i*N*N+j*N+k+1]     = room[i*N*N+j*N+k+1] + dCon[0];
				}
				if(j==N-1 || room[i*N*N+j*N+k+N] == -1){
					dCon[1] = 0;
				}else{
					dCon[1] = (room[i*N*N+j*N+k]    -room[i*N*N+j*N+k+N]) * coefficient;
					room[i*N*N+j*N+k] = room[i*N*N+j*N+k] - dCon[1];
					room[i*N*N+j*N+k+N]     = room[i*N*N+j*N+k+N] + dCon[1];
				}
				if(i==N-1 || room[i*N*N+j*N+k+(N*N)] == -1){
					dCon[2] = 0;
				}else{
					dCon[2] = (room[i*N*N+j*N+k]-room[i*N*N+j*N+k+(N*N)]) * coefficient;
					room[i*N*N+j*N+k] = room[i*N*N+j*N+k] - dCon[2];
					room[i*N*N+j*N+k+(N*N)] = room[i*N*N+j*N+k+(N*N)] + dCon[2];
				}
				if(k==0 || room[i*N*N+j*N+k-1] == -1){
					dCon[3] = 0;
				}else{
					dCon[3] = (room[i*N*N+j*N+k]    -room[i*N*N+j*N+k-1]) * coefficient;
					room[i*N*N+j*N+k] = room[i*N*N+j*N+k] - dCon[3];
					room[i*N*N+j*N+k-1]     = room[i*N*N+j*N+k-1] + dCon[3];
				}
				if(j==0 || room[i*N*N+j*N+k-N] == -1){
					dCon[4] = 0;
				}else{
					dCon[4] = (room[i*N*N+j*N+k]    -room[i*N*N+j*N+k-N]) * coefficient;
					room[i*N*N+j*N+k] = room[i*N*N+j*N+k] - dCon[4];
					room[i*N*N+j*N+k-N]     = room[i*N*N+j*N+k-N] + dCon[4];
				}
				if(i==0 || room[i*N*N+j*N+k-(N*N)] == -1){
					dCon[5] = 0;
				}else{
					dCon[5] = (room[i*N*N+j*N+k]-room[i*N*N+j*N+k-(N*N)]) * coefficient;
					room[i*N*N+j*N+k] = room[i*N*N+j*N+k] - dCon[5];
					room[i*N*N+j*N+k-(N*N)] = room[i*N*N+j*N+k-(N*N)] + dCon[5];
				}
				}
			}
		}
	}

	//after resetting the concentration values we then find the values of min and max
	//	in order to tell when the loop shall end
	conMin = room[0];
	conMax = room[0];

	for (i=0; i<N; i++) {
                for (j=0; j<N; j++){
                        for (k=0; k<N; k++){
				if (room[i*N*N+j*N+k] < conMin && room[i*N*N+j*N+k] != -1) {
                                        conMin = room[i*N*N+j*N+k];
                                }
                                if (room[i*N*N+j*N+k] > conMax && room[i*N*N+j*N+k] != -1) {
                                        conMax = room[i*N*N+j*N+k];
                                }	
			}
		}
	}

}
