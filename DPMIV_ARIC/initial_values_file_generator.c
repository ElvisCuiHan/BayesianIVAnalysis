#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>


int main(void)
{
	FILE  *temp;
	double init_b1=-9.862;
	// use posterior means from the IV-normal model;
	double init_a1[]={0.648};
	double init_a2[]={0.001,0.011,0.002,0.002,-0.004,-0.002,0.012};
	double init_b2[]={1.458,-0.300,-0.055,-0.599,1.629,-1.533,-1.668};
	double init_thetac[]={1.667, 0.008, 67.391, 19.543, 0.103};
	double init_v = 1.0;
	int init_ncl = 1;
	int n=12625, kG=1, kU=7; 
	int *init_ci, i;
	init_ci = (int *) malloc(n*sizeof(int)); // cluster indicator
	for (i=0; i<n; i++) init_ci[i]=0;

temp=fopen("temp_init_values.txt", "w");
			fprintf(temp, "%f\n", init_b1);
			for(i=0; i<kG; i++) fprintf(temp, "%f	", init_a1[i]); fprintf(temp,"\n");
			for(i=0; i<kU; i++) fprintf(temp, "%f	", init_a2[i]); fprintf(temp,"\n");
			for(i=0; i<kU; i++) fprintf(temp, "%f	", init_b2[i]); fprintf(temp,"\n");
			fprintf(temp, "%f\n", init_v);
			fprintf(temp, "%d\n", init_ncl); 
			for(i=0; i<5; i++) fprintf(temp, "%f	", init_thetac[i]); 
			fprintf(temp, "%d\n", n); 				
			for(i=0; i<n; i++) fprintf(temp, "%d	", init_ci[i]); fprintf(temp,"\n");
			fclose(temp);

  return 0;

}
