#include <stdio.h>
#include <mpi.h>

int main(int argc, char* argv[]) {
  int myrank, nproz;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproz);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if (myrank == 0) {
    printf("Size: %2d\n", nproz);
  }

  printf("Rank: %2d\n", myrank);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}
