#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>

#define N 10

int main(){
	double h, hsqrd, lRoom, tStep, dcoef, tAccum, cMax, cMin, urms;
	int rank, size, ierr, i, j, k, t0=0, t1=1, it;
	int partition;
	
	double *cube = malloc((N+2)*(N+2)*(N+2)*sizeof(double));
	double *mask = malloc((N+2)*(N+2)*(N+2)*sizeof(double));
    	double *cubeCopy = malloc((N)*(N)*(N)*sizeof(double));

	for(int i = 1; i <= N; i++){
                for(int j = 1; j <= N; j++){
                        for(int k =1; k <= N; k++){
                                mask[i,j,k] = 1.0;//mask[i*N*N+j*N+k] = 1.0;
                        }
                }
        }



    	lRoom = 5.0;
	urms = 250.0;
    	tStep = (lRoom/urms)/(double)N;
	h = lRoom/(double)N;
    	hsqrd = pow(h,2.0);

    	//cube[0] = 1e21;
    	/*dcoef = 0.175 * tStep / hsqrd;
	cMax = cube[1*N*N+1*N+1];
	cMin = 0;
	tAccum = 0;*/

	for(int i = 0; i < N+2; i++){
                for(int j = 0; j < N+2; j++){
                        for(int k =0; k < N+2; k++){
				cube[i,j,k] = 0.0;//cube[i*N*N+j*N+k] = 0.0;
			}
		}
	}
	cube[1*N*N+1*N+1] = 1e21;
	dcoef = 0.175 * tStep / hsqrd;
        cMax = cube[1,1,1];
        cMin = 0;
        tAccum = 0;

    	while(cMin/cMax <= 0.99){

		tAccum = tAccum + tStep;
		//printf("%f\n", cube[0]);
		
        	for(int i = 1; i <= N; i++){
            		for(int j = 1; j <= N; j++){
                		for(int k =1; k <= N; k++){

					//printf("%d %d %d\n", i, j, k);
					//double c_current_cell = C[i][j][k];
                    			/*cubeCopy[(i-1)*N*N+(j-1)*N+(k-1)] = cube[i*N*N+j*N+k] + (
                    			cube[(i+1)*N*N+j*N+k]*mask[(i+1)*N*N+j*N+k] + cube[(i-1)*N*N+j*N+k]*mask[(i-1)*N*N+j*N+k] +
                    			cube[i*N*N+(j+1)*N+k]*mask[i*N*N+(j+1)*N+k] + cube[i*N*N+(j-1)*N+k]*mask[i*N*N+(j-1)*N+k] +
                    			cube[i*N*N+j*N+(k+1)]*mask[i*N*N+j*N+(k+1)] + cube[i*N*N+j*N+(k-1)]*mask[i*N*N+j*N+(k-1)] -
                    			(double)(mask[(i+1)*N*N+j*N+k]+mask[(i-1)*N*N+j*N+k]+mask[i*N*N+(j+1)*N+k]+mask[i*N*N+(j-1)*N+k]+
						mask[i*N*N+j*N+(k+1)]+mask[i*N*N+j*N+(k-1)])*cube[i*N*N+j*N+k])*2*dcoef;
					*/

					cubeCopy[i-1,j-1,k-1] = cube[i,j,k] + (
                                        cube[(i+1),j,k]*mask[(i+1),j,k] + cube[(i-1),j,k]*mask[(i-1),j,k] +
                                        cube[i,(j+1),k]*mask[i,(j+1),k] + cube[i,(j-1),k]*mask[i,(j-1),k] +
                                        cube[i,j,(k+1)]*mask[i,j,(k+1)] + cube[i,j,(k-1)]*mask[i,j,(k-1)] -
                                        (double)(mask[(i+1),j,k]+mask[(i-1),j,k]+mask[i,(j+1),k]+mask[i,(j-1),k]+
                                                mask[i,j,(k+1)]+mask[i,j,(k-1)])*cube[i,j,k])*2*dcoef;

					
					/*if(cMin > c_current_cell){
		                            cMin = c_current_cell;
                		            mini = i;
                            		    minj = j;
                         		    mink = k;   
                        		}
                        		if(c_max < c_current_cell)
                            			cMax = c_current_cell;*/
                		}
            		}
        	}
		
		/*for (i=0; i<N; i++) {
	                for (j=0; j<N; j++){
        	                for (k=0; k<N; k++){
					if (cubeCopy[i*N*N+j*N+k] < cMin) {
                        	                cMin = cubeCopy[i*N*N+j*N+k];
                                	}
                                	if (cubeCopy[i*N*N+j*N+k] > cMax) {
                                        	cMax = cubeCopy[i*N*N+j*N+k];
                                	}	
				}
			}
		}*/
		
		cMax = cubeCopy[0];
		for (i=1;i<(N)*(N)*(N);i++)
			if(cubeCopy[i] > cMax)
				cMax = cubeCopy[i];
		cMin = cubeCopy[0];
		for (i=1;i<(N)*(N)*(N);i++)
                        if(cubeCopy[i] < cMin)
                                cMin = cubeCopy[i];


		//printf("%f\n", cMin);
		//printf("%f\n", cMax);
		//printf("%f\n", cMin/cMax);


		
		for(int i = 1; i <= N; i++){
                        for(int j = 1; j <= N; j++){
                                for(int k = 1; k <= N; k++){
                    			cube[i,j,k] = cubeCopy[(i-1),(j-1),(k-1)];
                    			//printf("%f\n", cube[i*N*N+j*N+k]);
                		}
            		}
        	}


		double tot = 0;
	        for (i=0; i<N+2; i++) {
        	        for (j=0; j<N+2; j++){
                	        for (k=0; k<N+2; k++){
                        	        tot = tot + cube[i,j,k];
                        	}
                	}
        	}
		printf("%f\n",tot);

    	}
	
	double tot = 0;
	for (i=0; i<N+2; i++) {
                for (j=0; j<N+2; j++){
                        for (k=0; k<N+2; k++){
                                tot = tot + cube[i,j,k];
                        }
                }
        }

    	printf("%f\n",tot);
	printf("%f\n",tAccum);
    	free(cube);
	free(cubeCopy);
	free(mask);
}
