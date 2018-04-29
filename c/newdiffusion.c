#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>
#include <sys/time.h>
#include <time.h>

double factor = 1.0e-6;

// C Code checked on 10/10/17 for mass consistency

//Code given from email that will calculate wall time
double walltime_() 
{
    struct timeval tp;
    int rtn;
    double seconds;

    rtn=gettimeofday(&tp, NULL);

    seconds = tp.tv_sec + factor * tp.tv_usec;

    return  seconds ; 
}
void diffuse(double*** C, double*** Ccopy, int*** mask, int M, int partition);
int main(int argc, char** argv)
{
      double cpu;
      const int M;
      printf("What is the size of the box?\n");
      //Variable declarations
      scanf("%d", &M);
      if(M <= 0)
      {
            printf("Invalid number");
            exit(0);
      }
      int partition;
      printf("Is there a partition? (0 for no, 1 for yes)\n");
      scanf("%d", &partition);
      if(partition != 1 && partition != 0)
      {
            printf("Partition will be turned off.\n");
            partition = 0;
      }
      int i,j,k;
      //Makes the array into cube form
      double ***C = malloc((M+2)*sizeof(double**));
      
      for(i = 0; i < M+2; i++)
      {
            C[i] = malloc((M+2)*sizeof(double*));
            for(j = 0; j < M+2; j++)
            {
                  C[i][j] = malloc((M+2)*sizeof(double));
            }
      }
      //Zeroes out the array
      for(i = 0; i < M+2; i++)
      {
            for(j = 0; j < M+2; j++)
            {
                  for(k = 0; k < M+2; k++)
                  {
                        C[i][j][k] = 0.0;
                  }
            }
      }
      double ***Ccopy = malloc((M+2)*sizeof(double**));
      for(i = 0; i < M+2; i++)
      {
            Ccopy[i] = malloc((M+2)*sizeof(double*));
            for(j = 0; j < M+2; j++)
            {
                  Ccopy[i][j] = malloc((M+2)*sizeof(double));
            }
      }
      //Zeroes out the array
      for(i = 0; i < M+2; i++)
      {
            for(j = 0; j < M+2; j++)
            {
                  for(k = 0; k < M+2; k++)
                  {
                        Ccopy[i][j][k] = 0.0;
                  }
            }
      }
      
      int ***mask = malloc((M+2)*sizeof(int**));
      for(i = 0; i < M+2; i++)
      {
            mask[i] = malloc((M+2)*sizeof(int*));
            for(j = 0; j < M+2; j++)
            {
                  mask[i][j] = malloc((M+2)*sizeof(int));
            }
      }
      //Zeroes out the array
      for(i = 0; i < M+2; i++)
      {
            for(j = 0; j < M+2; j++)
            {
                  for(k = 0; k < M+2; k++)
                  {
                        mask[i][j][k] = 0;
                  }
            }
      }
      //Set the inside of cube to 1
      for(i = 1; i < M+1; i++)
      {
            for(j = 1; j < M+1; j++)
            {
                  for(k = 1; k < M+1; k++)
                  {
                        mask[i][j][k] = 1;
                  }
            }
      }
      
      //Starts the counting
      cpu = walltime_();
      printf("Beginning Box Simulation...\n");
      diffuse(C, Ccopy, mask, M, partition);  //Calls the method diffuse
      free(C);     //Empties C to save space
      free(Ccopy);
      free(mask);
      cpu = walltime_() - cpu;
      printf("CPU time = %f\n", cpu);
}
//Method to go through the array and diffuse the box
void diffuse(double*** C, double*** Ccopy, int*** mask, int M, int partition)
{
      int partsize;
      if(partition == 1)
      {
            partsize = floor(M / 2);
      }
      //More variables
      C[1][1][1] = 1.0 * pow(10,21);  //Makes first cell 1.0e21
      double room = 5.0;
      double diff = 0.175;
      double urms = 250.0;  //Given constant for g/mol of gas
      double tstep = (double) (room / urms) / M;
      double height = (double) room / M;
      double tacc = 0.0;
      double ratio = 0.0;
      double sum = 0.0;
      double dC = diff * tstep / (height * height);
      double change;
      if(partition == 1)
      {
            for(int i = 0; i < M; i++)
            {
                  for(int j = 0; j < M; j++)
                  {
                        for(int k = 0; k < M; k++)
                        {
                              if((i == partsize-1) && (j >= partsize-1))
                              {
                                    C[i][j][k] = -1.0;
                              }
                        }
                  }
            }      
      }
      printf("Before the loop\n");
      //Loop that checks if the boxes are not all equal
      do
      {
            for(int i = 1; i < M+1; i++)
            {
                  for(int j = 1; j < M+1; j++)
                  {
                        for(int k = 1; k < M+1; k++)
                        {
                              //Compare this to equation in book.  I also added
                              //the diffusion coefficient which you were missing.
                              Ccopy[i][j][k] = C[i][j][k] + (
                              C[i+1][j][k]*mask[i+1][j][k] + C[i-1][j][k]*mask[i-1][j][k] + 
                              C[i][j+1][k]*mask[i][j+1][k] + C[i][j-1][k]*mask[i][j-1][k] + 
                              C[i][j][k+1]*mask[i][j][k+1] + C[i][j][k-1]*mask[i][j][k-1] - 
                              (double)(mask[i+1][j][k]+mask[i-1][j][k]+mask[i][j-1][k]+mask[i][j+1][k] +
                              mask[i][j][k-1]+mask[i][j][k+1])*C[i][j][k]) * 2 * dC;
                        }
                  }
            }
            //printf("After Jacobi\n");
            tacc = tacc + tstep;
            //Makes sure there is mass consistency
            double maxc = C[1][1][1];
            double minc = C[1][1][1];
            sum = 0.0;
            for(int i = 0; i < M+2; i++)
            {
                  for(int j = 0; j < M+2; j++)
                  {
                        for(int k = 0; k < M+2; k++)
                        {
                              if(C[i][j][k] != -1.0)
                              {
                                    maxc = fmax(C[i][j][k], maxc);
                                    minc = fmin(C[i][j][k], minc);
                                    sum += Ccopy[i][j][k];
                              }
                        }
                  }
            }
            ratio = minc / maxc;      //Sees if min and max are equal in order to end while loop
            //Prints the different variable types in loop
              //printf("%f %f %f\n", tacc, ratio, C[0][0][0]);
              //printf("%f\n", C[M-1][M-1][M-1]);
              //printf("%f\n", sum);
      } while(ratio <= 0.99);
      printf("Total sum is %f.\n", sum);
      printf("Box completed in %f seconds.\n", tacc);
}
