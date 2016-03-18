#include <stdio.h>
#include <mpi.h>

int main(int argc, char** argv) {
  int nranks, rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    float pi = 22.0 / 7.0;

    for (int i = 1; i < nranks; i++) {
      printf("Send to rank %2d\n", i);
      MPI_Send(&pi, 1, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
    }
  } else {
    float pi;

    MPI_Recv(&pi, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    printf("Rank %2d -> pi = %f\n", rank, pi);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}
