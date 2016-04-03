#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

void printvec(double *x, int n) {
  printf("[");
  for (int i = 0; i < n; i++) {
    if (i > 0) {
      printf(", ");
    }

    printf("%f", x[i]);
  }
  printf("]\n");
}

void printmat(double *A, int n, int m) {
  for (int i = 0; i < n; i++) {
    printvec(A + i * m, m);
  }
}

int main(int argc, char **argv) {
  int rank, nranks;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);

  srand(5);

  int kstop = 20;
  double omega = 1.5;
  int n = 4000;
  int N = n / nranks;
  double *A = malloc(N * n * sizeof(A[0]));
  double *b = malloc(n * sizeof(b[0]));

  if (rank == 0) {
    double *A2 = malloc(n * n * sizeof(A2[0]));

    A2[0] = 4.0;
    A2[1] = -1.0;

    for (int i = 1; i < n - 1; i++) {
      A2[i * n + i - 1] = -1.0;
      A2[i * n + i] = 4.0;
      A2[i * n + i + 1] = -1.0;
    }

    A2[n * n - 2] = -1.0;
    A2[n * n - 1] = 4.0;

    for (int i = 0; i < n; i++) {
      b[i] = (double)rand() / RAND_MAX;
    }

    for (int j = 0; j < n / nranks; j++) {
      memcpy(A + j * n, A2 + j * nranks * n, n * sizeof(double));
    }

    for (int i = 1; i < nranks; i++) {
      for (int j = 0; j < n / nranks; j++) {
        MPI_Send(A2 + (i + j * nranks) * n, n, MPI_DOUBLE, i, 0,
                 MPI_COMM_WORLD);
      }
    }
  } else {
    for (int i = 0; i < N; i++) {
      MPI_Recv(A + i * n, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
    }
  }

  MPI_Bcast(b, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  double *alpha = malloc(n * sizeof(alpha[0]));
  double *c = malloc(N * n * sizeof(c[0]));
  double *x = malloc(n * sizeof(x[0]));

  int col = 0;
  for (int i = 0; i < n; i++) {
    if (i % nranks == rank) {
      alpha[i] = omega / A[col * n + i];
      col += 1;
    }

    x[i] = 0.0;
  }

  for (int j = 0; j < N; j++) {
    int col = j * nranks + rank;

    for (int i = 0; i < n; i++) {
      if (i == col) {
        c[j * n + i] = (1.0 - omega) * A[j * n + i] / omega;
      } else {
        c[j * n + i] = -A[j * n + i];
      }
    }
  }

  for (int k = 0; k < kstop; k++) {
    for (int i = 0; i < n; i++) {
      double a = 0.0;

      for (int j = 0; j < N; j++) {
        int row = j * nranks + rank;

        a += c[j * n + i] * x[row];
      }

      if (i % nranks == rank) {
        double *as = malloc(nranks * sizeof(as[0]));

        MPI_Gather(&a, 1, MPI_DOUBLE, as, 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);

        double s = 0.0;

        for (int k = 0; k < nranks; k++) {
          s += as[k];
        }

        free(as);

        x[i] = alpha[i] * (s + b[i]);
      } else {
        MPI_Gather(&a, 1, MPI_DOUBLE, NULL, nranks, MPI_DOUBLE, i % nranks,
                   MPI_COMM_WORLD);
      }
    }
  }

  if (rank == 0) {
    double buf;

    for (int i = 0; i < n; i++) {
      if (i % nranks != 0) {
        MPI_Recv(x + i, 1, MPI_DOUBLE, i % nranks, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
      }
    }
  } else {
    for (int i = 0; i < n; i++) {
      if (i % nranks == rank) {
        MPI_Send(x + i, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      }
    }
  }

  if (rank == 0) {
    printvec(x, n);
  }

  free(x);
  free(c);
  free(alpha);

  free(b);
  free(A);

  MPI_Finalize();

  return 0;
}
