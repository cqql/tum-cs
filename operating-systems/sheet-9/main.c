#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// generate random matrix size NxN
float ** create_random_matrix(int n) {
	// allocate columns
	float **matrix =  malloc(n*sizeof(float*));
	int i,j;

	for(i = 0; i < n; i++) {
		// allocate rows
		matrix[i] = malloc(n*sizeof(float));
		for(j = 0; j < n; j++) {
			matrix[i][j] = (float) rand() / (float) RAND_MAX;
		}
	}

	return matrix;
}

// free the previously allocated matrix
void free_matrix(float** matrix, int n) {
	int i;

	for(i = 0; i < n; i++) {
		// free rows
		free( (void*) matrix[i] );
	}

	free( (void*) matrix );
}

/* __________________________________
/ Code for timing:                   \
| use tic() to start                  |
| measurement and                     |
| double elapsed_seconds=toc()           |
\ to end measurement                 /
 ------------------------------------
   \
    \
        .--.
       |o_o |
       |:_/ |
      //   \ \
     (|     | )
    /'\_   _/`\
    \___)=(___/
*/

// current cpu-time stamp value
clock_t cpu_start_time;

struct timespec usr_ts1;

// start cpu-time measurement
void tic() {
	clock_gettime(CLOCK_REALTIME, &usr_ts1);
}

// end cpu-time measurement (seconds)
double toc() {
	double usr_elapsed;
        struct timespec usr_ts2;

  	clock_gettime(CLOCK_REALTIME, &usr_ts2);
  	usr_elapsed = usr_ts2.tv_sec + usr_ts2.tv_nsec/1E9;
  	usr_elapsed -= usr_ts1.tv_sec + usr_ts1.tv_nsec/1E9;

	return usr_elapsed;
}

int main() {
	// set a new seed for the random number generator
	// at each start (system time)
	srand(time(NULL));

	//////////////////////////
	//
	// assignment code  starts here :)
	//
	//////////////////////////

	int n = 4096;
  int N = 20;
	// create random matrix and measure the time elapsed for
	// creation
	tic();
	float** matrix = create_random_matrix(n);
	double elapsed_time = toc();
	printf("Elapsed time creating a random matrix: %f seconds\n",elapsed_time);

  double sum = 0;
  elapsed_time = 0;
  for (int m = 0; m < N; m++) {
    tic();

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        sum += matrix[i][j];
      }
    }

    elapsed_time += toc();
  }
  sum /= N;
  elapsed_time /= N;
  printf("row-wise %f seconds -> %f\n", elapsed_time, sum);

  sum = 0;
  elapsed_time = 0;
  for (int m = 0; m < N; m++) {
    tic();

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        sum += matrix[j][i];
      }
    }

    elapsed_time += toc();
  }
  sum /= N;
  elapsed_time /= N;
  printf("column-wise %f seconds -> %f\n", elapsed_time, sum);

  int n2 = n * n;
  int* indices = malloc(n2 * sizeof(int));

  for (int i = 0; i < n2; i++) {
    indices[i] = (rand() / RAND_MAX) * n;
  };

  sum = 0;
  elapsed_time = 0;
  for (int m = 0; m < N; m++) {
    tic();

    for (int i = 0; i < n2; i++) {
      sum += matrix[indices[i]][indices[n2 - i - 1]];
    }

    elapsed_time += toc();
  }
  sum /= N;
  elapsed_time /= N;
  printf("random %f seconds -> %f\n", elapsed_time, sum);

  free(indices);

	free_matrix(matrix, n);

	return 0;
}
