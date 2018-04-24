#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>

#define N 10

int main(){
	double deltas, lRoom, tStep, dcoef;
	int rank, size, ierr, i, j, k, t0=0, t1=1, it;
	int partition;
	
	double *cube = malloc((N+1)*(N+1)*(N+1)*sizeof(double));
    	double *cubeCopy = malloc(N*N*N*sizeof(double));
    	lRoom = 5;
    	tStep = (lRoom/250.0)/(double)N;
    	deltas = lRoom/(double)N;
    	deltas = pow(deltas,2.0);

    	//cube[0] = 1e21;
    	dcoef = 0.175 * tStep / deltas;

	for(int i = 0; i < N+1; i++){
                for(int j = 0; j < N+1; j++){
                        for(int k =0; k < N+1; k++){
				cube[i*N*N+j*N+k] = 0.0;
			}
		}
	}
	cube[N*N+N+1] = 1e21;


    	for(int it =1; it < 1000; it++){

		//printf("%f\n", cube[0]);

        	for(int i = 1; i < N; i++){
            		for(int j = 1; j < N; j++){
                		for(int k =1; k < N; k++){

					//printf("%d %d %d\n", i, j, k);

                    			cubeCopy[i*N*N+j*N+k] = cube[i*N*N+j*N+k] + (
                    			cube[(i+1)*N*N+j*N+k] + cube[(i-1)*N*N+j*N+k] +
                    			cube[i*N*N+(j+1)*N+k] + cube[i*N*N+(j-1)*N+k] +
                    			cube[i*N*N+j*N+(k+1)] + cube[i*N*N+j*N+(k-1)] -
                    			6.0*(cube[i*N*N+j*N+k])) * dcoef / deltas;
                		}
            		}
        	}
		for(int i = 1; i < N; i++){
                        for(int j = 1; j < N; j++){
                                for(int k = 1; k < N; k++){
                    			cube[i*N*N+j*N+k] = cubeCopy[i*N*N+j*N+k];
                    			printf("%f\n", cube[i*N*N+j*N+k]);
                		}
            		}
        	}
    	}
	
	double tot = 0;
	for (i=0; i<N+1; i++) {
                for (j=0; j<N+1; j++){
                        for (k=0; k<N+1; k++){
                                tot = tot + cube[i*N*N+j*N+k];
                        }
                }
        }

    	printf("%f\n",tot);
    	free(cube);
}
