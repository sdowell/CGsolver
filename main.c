/* UCSB CS240A, Winter Quarter 2014
 * Main and supporting functions for the Conjugate Gradient Solver on a 5-point stencil
 *
 * NAMES:
 * PERMS:
 * DATE:
 */
#include "mpi.h"
#include "hw2harness.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double* load_vec( char* filename, int* k );
void save_vec( int k, double* x );
double ddot(double *v, double *w, int n);
double* daxpy(double alpha, double *v, double beta, double *w, int n);
double* matvec(double *v, int k, int rank, int size);
double* cgsolve(double *b, int k, int rank, int size,int maxiters, double relres);

int main( int argc, char* argv[] ) {
	int writeOutX = 0;
	int n, k;
	int maxiterations = 100000;
	int niters=0;
 	double norm;
	double* b;
	double* x;
	double time;
	double t1, t2;
	
	MPI_Init( &argc, &argv );
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	// Read command line args.
	// 1st case runs model problem, 2nd Case allows you to specify your own b vector
	if ( argc == 3 ) {
		k = atoi( argv[1] );
		n = k*k;
		// each processor calls cs240_getB to build its own part of the b vector!
	} else if  ( !strcmp( argv[1], "-i" ) && argc == 4 ) {
		b = load_vec( argv[2], &k );
	} else {
		printf( "\nCGSOLVE Usage: \n\t"
			"Model Problem:\tmpirun -np [number_procs] cgsolve [k] [output_1=y_0=n]\n\t"
			"Custom Input:\tmpirun -np [number_procs] cgsolve -i [input_filename] [output_1=y_0=n]\n\n");
		exit(0);
	}
	writeOutX = atoi( argv[argc-1] ); // Write X to file if true, do not write if unspecified.

	
	// Start Timer
	t1 = MPI_Wtime();
	
	// CG Solve here!
        double x_initial[n];
	x = x_initial;
	double *r;
	double *rnew;
	double *d;
	double alpha;
	double beta;
	
	// Test to check values of B
/*	if(rank == 0){
		int i;
		for(i = 0; i < n; i++){
			printf("b[%d]: %f\n", i, cs240_getB(i,n));
		}
	}	
*/	
	// Test to exercise distributing ddot across multiple processors
	int myN = n / size;
	//printf("myN: %d = %d / %d\n", myN, n, size);
	//double *tmp = malloc(sizeof(double) * myN);
	//double tmp[myN];
	int i;
/*	for(i = 0; i < myN; i++){
		tmp[i] = cs240_getB(i + (rank * myN), n);
		printf("Proc %d: tmp[%d] = %f\n", rank, i + (rank * myN), tmp[i]);
	}
	double dd = ddot(tmp, tmp, myN);
	printf("Ddot from process %d: %lf\n", rank, dd);
*/	
/*	printf("-------------------Testing matvec-----------------------\n");
	for(i = 0; i < myN; i++){
		tmp[i] = cs240_getB(i + (rank * myN), n);
	}	
	double *v = matvec(tmp, k, rank, size);
	int ind = rank * myN;
	for(i = 0; i < myN; i++){
		printf("Proc %d: v[%d] = %f\n", rank, ind + i, v[i]);
	}
*/
/*	r = b;
	d = r;
	while(){
		alpha = ddot(r, r) / matvec();
		x = x + alpha 
	}
*/
	b = malloc(sizeof(double) * myN);
	for(i = 0; i < myN; i++){
		b[i] = cs240_getB(i + (rank * myN), n);
	}
	x = cgsolve(b, k, rank, size, maxiterations, 1.0);	
	//int correct = cs240_verify(x, k, 1.0);
	//printf("P%d Correct: %d\n",rank, correct);
 	// End Timer
	t2 = MPI_Wtime();
	int correct = cs240_verify(x, k, 1.0);
	printf("P%d Correct: %d\n", rank, correct);	
	if ( writeOutX ) {
		save_vec( k, x );
	}
		
	// Output
	printf( "Problem size (k): %d\n",k);
	if(niters>0){
          printf( "Norm of the residual after %d iterations: %lf\n",niters,norm);
        }
	printf( "Elapsed time during CGSOLVE: %lf\n", t2-t1);
	
        // Deallocate 
        if(niters > 0){
	  free(b);
	}
        if(niters > 0){
          free(x);
	}
	
	MPI_Finalize();
	
	return 0;
}
double* cgsolve(double *b, int k, int rank, int size, int maxiters, double relres){
//	if(rank == 0)
//		printf("Entering matvec\n");
	relres = 1;
	int niters = 0;
	int n = k * k;
	int p = size;
	int myN = n / p;
	double *x = malloc(sizeof(double) * myN);
	double *r = malloc(sizeof(double) * myN);
//	if(rank == 0)
//		printf("memcpy r,b\n");
	memcpy(r, b, myN * sizeof(double));
	double rtr = ddot(r, r, myN);
	double *d = malloc(sizeof(double) * myN);
//	if(rank == 0)
//		printf("memcpy d,r\n");
	memcpy(d, r, myN * sizeof(double));
	double normb = sqrt(ddot(b, b, myN));
	//double *Ad;
	double alpha, beta;
	double rtrold;
//	if(rank == 0)
//		printf("Beginning matvec loop\n");
	while(relres > .000001 && niters < maxiters){
		niters++;
		double *Ad = matvec(d, k, rank, size);
		alpha = rtr / ddot(d, Ad, myN);
		daxpy(1.0, x, alpha, d, myN);
		daxpy(1.0, r, -1.0 * alpha, Ad, myN);
		rtrold = rtr;
		rtr = ddot(r, r, myN);
		beta = rtr / rtrold;
		daxpy(beta, d, 1, r, myN);
		relres = sqrt(rtr) / normb;
		free(Ad);
	}
//	return x;
	double *xbuf = malloc(sizeof(double) * n);
	MPI_Allgather(x, myN, MPI_DOUBLE, xbuf, myN, MPI_DOUBLE, MPI_COMM_WORLD);
/*	if(rank == 0){
		int i;
		for(i = 0; i < n; i++)
			printf("x[%d] = %f\n", i, x[i]);
	}
*/	free(x);
	free(r);
	free(d);
	return xbuf;
}
double ddot(double *v, double *w, int n){
	double localsum = 0;
	double globalsum = 0;
	
	int i;
	for(i = 0; i < n; i++){
		localsum += v[i] * w[i];
		//printf("Proc %d: %f += %f * %f\n",i, localsum, v[i], w[i]);
	}
	
	MPI_Allreduce(&localsum, &globalsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return globalsum;
}

double* daxpy(double alpha, double *v, double beta, double *w, int n){
	int i;
	for(i = 0; i < n; i++){
		v[i] = (alpha * v[i]) + (beta * w[i]);
	}
}
double* matvec(double *w, int k, int rank, int size){
	int n = k * k;
	int p = size;
	double *v = malloc(sizeof(double) * n / p);
	MPI_Status status;
	double *ww = malloc(sizeof(double) * (n / p + 2 * k));
	memcpy(&ww[k], w, (n / p) * sizeof(double));
	int i, j;
/*	for(i = 0; i < n / p; i++){
		if(rank == 0)
			printf("Proc %d: w[%d] = %f\n", rank, i, w[i]);
	}
	for(i = 0; i < n / p + 2 * k; i++){
		if(rank == 0){
			printf("Proc %d: ww[%d] = %f\n", rank, i, ww[i]);
		}
	}
*/	// TODO: construct boundary conditions for w from adjacent processors
	if(rank % 2 == 1){
		// TODO: write then read from previous processor
		MPI_Send(&w[0], k, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
		MPI_Recv(&ww[0], k, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
		// TODO: write then read from next processor
		if(rank != size - 1){
			MPI_Send(&w[n / p - k], k, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
			MPI_Recv(&ww[n / p + k], k, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status);
		}
	}else{
		// TODO: read then write from next processor
		if(rank != size - 1){
			MPI_Recv(&ww[n / p + k], k, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status);
			MPI_Send(&w[n / p - k], k, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
		}
		// TODO: read then write from previous processor
		if(rank != 0){
			MPI_Recv(&ww[0], k, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
			MPI_Send(&w[0], k, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);	
		}
	}
	//int i, j;
/*	for(i = 0; i < n / p + 2 * k; i++){
		if(rank == 0)
			printf("Proc %d: ww[%d] = %f\n", rank, i, ww[i]);
	}
*/	double *wptr = &ww[k];
	int start = k * rank / size;
	int end = start + k / size;
	int globalindex, localindex;
	for(i = start; i < end; i++){
		for(j = 0; j < k; j++){
			globalindex = (i) * k + j;
			localindex = globalindex - start * k;
			v[localindex] = 4 * wptr[localindex];
			if(i != 0)
				v[localindex] -= wptr[localindex - k];
			if(j != 0)
				v[localindex] -= wptr[localindex - 1];
			if(j != k - 1)
				v[localindex] -= wptr[localindex + 1];
			if(i != k - 1)
				v[localindex] -= wptr[localindex + k];
		}
	}
	free(ww);
	return v;
}

int getA(int index, int n){

}

/*
 * Supporting Functions
 *
 */

// Load Function
// NOTE: does not distribute data across processors
double* load_vec( char* filename, int* k ) {
	FILE* iFile = fopen(filename, "r");
	int nScan;
	int nTotal = 0;
	int n;
	
	if ( iFile == NULL ) {
		printf("Error reading file.\n");
		exit(0);
	}
	
	nScan = fscanf( iFile, "k=%d\n", k );
	if ( nScan != 1 ) {
		printf("Error reading dimensions.\n");
		exit(0);
	}
	
	n = (*k)*(*k);
	double* vec = (double *)malloc( n * sizeof(double) );
	
	do {
		nScan = fscanf( iFile, "%lf", &vec[nTotal++] );
	} while ( nScan >= 0 );
	
	if ( nTotal != n+1 ) {
		printf("Incorrect number of values scanned n=%d, nTotal=%d.\n",n,nTotal);
		exit(0);
	}
	
	return vec;
}

// Save a vector to a file.
void save_vec( int k, double* x ) { 
	FILE* oFile;
	int i;
	oFile = fopen("xApprox.txt","w");
	
	fprintf( oFile, "k=%d\n", k );
	
	for (i = 0; i < k*k; i++) { 
    	fprintf( oFile, "%lf\n", x[i]);
 	} 

	fclose( oFile );
}


