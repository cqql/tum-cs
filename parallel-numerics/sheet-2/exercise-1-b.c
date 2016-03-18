#include <stdio.h>
#include <math.h>

#include <mpi.h>

double f(double x) {
  return -3 * pow(x, 2) + 3;
}

double integrate(double (*f)(double), double a, double b, int n) {
  double h = (b - a) / n;
  double sum = 0.0;

  for (int i = 0; i < n; i++) {
    sum += (f(a + i * h) + f(a + (i + 1) * h)) * h / 2;
  }

  return sum;
}

int main(int argc, char** argv) {
  int nranks, rank;
  double a = -1;
  double b = 1;
  int n;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc == 1) {
    n = 2;
  } else {
    sscanf(argv[1], "%d", &n);
  }

  double h = (b - a) / nranks;

  /* Local settings */
  a = a + rank * h;
  b = a + h;

  double sum = integrate(f, a, b, n);

  if (rank == 0) {
    for (int i = 1; i < nranks; i++) {
      double buf;
      MPI_Recv(&buf, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      sum += buf;
    }

    printf("Integral = %.10f\n", sum);
  } else {
    MPI_Send(&sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }

  MPI_Finalize();
}
