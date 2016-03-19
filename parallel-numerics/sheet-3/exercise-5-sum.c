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

  float *a, *b, *r;
  int num = n / nranks;

  if (rank == 0) {
    num += n - (n / nranks) * nranks;
  }

  if (rank == 0) {
    a = (float *)malloc(n * sizeof(float));
    b = (float *)malloc(n * sizeof(float));
    r = (float *)malloc(n * sizeof(float));

    for (int i = 0; i < n; i++) {
      a[i] = (float)rand() / (float)RAND_MAX;
      b[i] = (float)rand() / (float)RAND_MAX;
    }

    int offset = num;
    for (int i = 1; i < nranks; i++) {
      MPI_Send(&a[offset], n / nranks, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
      MPI_Send(&b[offset], n / nranks, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
      offset += n / nranks;
    }
  } else {
    a = (float *)malloc(num * sizeof(float));
    b = (float *)malloc(num * sizeof(float));
    r = (float *)malloc(num * sizeof(float));

    MPI_Recv(a, num, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(b, num, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  double start = MPI_Wtime();

  for (int i = 0; i < num; i++) {
    r[i] = a[i] + b[i];
  }

  double end = MPI_Wtime();

  free(a);
  free(b);

  if (rank == 0) {
    int offset = num;
    for (int i = 1; i < nranks; i++) {
      MPI_Recv(&r[offset], n / nranks, MPI_FLOAT, i, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      offset += n / nranks;
    }

    float sum = 0.0;
    for (int i = 0; i < n; i++) {
      sum += r[i];
    }
    printf("Sum: %f\n", sum);
  } else {
    MPI_Send(r, num, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
  }

  free(r);

  if (rank == 0) {
    printf("Time: %fs\n", end - start);
  }

  MPI_Finalize();

  return 0;
}
