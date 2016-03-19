#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

int main(int argc, char **argv) {
  int nranks, rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int n = 1000;

  if (argc > 1) {
    sscanf(argv[1], "%d", &n);
  }

  float *a, *b;
  int num = n / nranks;

  if (n % nranks == 0) {
    num = n / nranks;
  } else {
    if (rank == nranks - 1) {
      num = n - (n / nranks + 1) * (nranks - 1);
    } else {
      num = n / nranks + 1;
    }
  }

  if (rank == 0) {
    a = (float *)malloc(n * sizeof(float));
    b = (float *)malloc(n * sizeof(float));

    srand(1);
    for (int i = 0; i < n; i++) {
      a[i] = (float)rand() / (float)RAND_MAX;
      b[i] = (float)rand() / (float)RAND_MAX;
    }
  } else {
    a = (float *)malloc(num * sizeof(float));
    b = (float *)malloc(num * sizeof(float));
  }

  MPI_Scatter(a, num, MPI_FLOAT, a, num, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Scatter(b, num, MPI_FLOAT, b, num, MPI_FLOAT, 0, MPI_COMM_WORLD);

  float z = 0.0;

  double start = MPI_Wtime();

  for (int i = 0; i < num; i++) {
    z += a[i] * b[i];
  }

  double end = MPI_Wtime();

  free(a);
  free(b);

  float total;
  MPI_Reduce(&z, &total, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    printf("Sum: %f\n", total);
    printf("Time: %fs\n", end - start);
  }

  MPI_Finalize();

  return 0;
}
