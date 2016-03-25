#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

void printvec(float *x, int n) {
  printf("[");
  for (int i = 0; i < n; i++) {
    if (i > 0) {
      printf(", ");
    }

    printf("%f", x[i]);
  }
  printf("]\n");
}

float *solve(int n, float *a, float *b, float *c, float *d) {}

int main(int argc, char **argv) {
  int n = 5;

  float A[] = {0.0, -1.0, -1.0, -1.0, -1.0};
  float B[] = {2.0, 2.0, 2.0, 2.0, 2.0};
  float C[] = {-1.0, -1.0, -1.0, -1.0, 0.0};
  float D[] = {2.5, 1.5, 2.0, 0.1, 5.0};

  int rank, nranks;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);

  int *sendcounts = NULL;
  int *displs = NULL;
  int count = n / nranks;

  if (rank < n % nranks) {
    count += 1;
  }

  float *a = malloc(count * sizeof(float));
  float *b = malloc(count * sizeof(float));
  float *c = malloc(count * sizeof(float));
  float *d = malloc(count * sizeof(float));

  if (rank == 0) {
    sendcounts = malloc(nranks * sizeof(int));
    displs = malloc(nranks * sizeof(int));

    for (int i = 0; i < nranks; i++) {
      sendcounts[i] = n / nranks;

      if (i < n % nranks) {
        sendcounts[i] += 1;
      }

      if (i == 0) {
        displs[i] = 0;
      } else {
        displs[i] = displs[i - 1] + sendcounts[i - 1];
      }
    }
  }

  MPI_Scatterv(A, sendcounts, displs, MPI_FLOAT, a, count, MPI_FLOAT, 0,
               MPI_COMM_WORLD);
  MPI_Scatterv(B, sendcounts, displs, MPI_FLOAT, b, count, MPI_FLOAT, 0,
               MPI_COMM_WORLD);
  MPI_Scatterv(C, sendcounts, displs, MPI_FLOAT, c, count, MPI_FLOAT, 0,
               MPI_COMM_WORLD);
  MPI_Scatterv(D, sendcounts, displs, MPI_FLOAT, d, count, MPI_FLOAT, 0,
               MPI_COMM_WORLD);

  /* Eliminate subdiagonal */
  for (int i = 1; i < count; i++) {
    float t = -a[i] / b[i - 1];

    if (rank == 0) {
      a[i] = 0;
    } else {
      a[i] = t * a[i - 1];
    }

    b[i] += t * c[i - 1];
    d[i] += t * d[i - 1];
  }

  /* Eliminate superdiagonal */
  for (int i = count - 3; i >= 0; i--) {
    float t = -c[i] / b[i + 1];
    a[i] += t * a[i + 1];
    c[i] = t * c[i + 1];
    d[i] += t * d[i + 1];
  }

  /* Eliminate c[i] of last row */
  float la, lb, lc, ld;

  if (nranks > 1) {
    if (rank == 0) {
      MPI_Recv(&la, 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      MPI_Recv(&lb, 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      MPI_Recv(&lc, 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      MPI_Recv(&ld, 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
    } else if (rank == nranks - 1) {
      MPI_Send(&a[0], 1, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD);
      MPI_Send(&b[0], 1, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD);
      MPI_Send(&c[0], 1, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD);
      MPI_Send(&d[0], 1, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD);
    } else if (rank % 2 == 0) {
      MPI_Recv(&la, 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      MPI_Recv(&lb, 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      MPI_Recv(&lc, 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      MPI_Recv(&ld, 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      MPI_Send(&a[0], 1, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD);
      MPI_Send(&b[0], 1, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD);
      MPI_Send(&c[0], 1, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD);
      MPI_Send(&d[0], 1, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD);
    } else {
      MPI_Send(&a[0], 1, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD);
      MPI_Send(&b[0], 1, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD);
      MPI_Send(&c[0], 1, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD);
      MPI_Send(&d[0], 1, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD);
      MPI_Recv(&la, 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      MPI_Recv(&lb, 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      MPI_Recv(&lc, 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      MPI_Recv(&ld, 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
    }
  }

  if (rank < nranks - 1) {
    float t = -c[count - 1] / lb;
    b[count - 1] += t * la;
    c[count - 1] = t * lc;
    d[count - 1] += t * ld;
  }

  /* Eliminate elements below main diagonal */
  float ub, uc, ud;

  if (nranks > 1) {
    if (rank == 0) {
      MPI_Send(&b[count - 1], 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD);
      MPI_Send(&c[count - 1], 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD);
      MPI_Send(&d[count - 1], 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD);
    } else {
      MPI_Recv(&ub, 1, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      MPI_Recv(&uc, 1, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      MPI_Recv(&ud, 1, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);

      float t = -a[count - 1] / ub;
      a[count - 1] = 0.0;
      b[count - 1] += t * uc;
      d[count - 1] += t * ud;

      if (rank < nranks - 1) {
        MPI_Send(&b[count - 1], 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD);
        MPI_Send(&c[count - 1], 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD);
        MPI_Send(&d[count - 1], 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD);
      }

      for (int i = 0; i < count - 1; i++) {
        t = -a[i] / ub;
        a[i] = 0.0;
        c[i] += t * uc;
        d[i] += t * ud;
      }
    }
  }

  /* Eliminate upper triagonal elements */
  if (nranks > 1) {
    if (rank == nranks - 1) {
      MPI_Send(&b[count - 1], 1, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD);
      MPI_Send(&d[count - 1], 1, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD);
    } else {
      MPI_Recv(&lb, 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      MPI_Recv(&ld, 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);

      float t = -c[count - 1] / lb;
      d[count - 1] += t * ld;

      if (rank > 0) {
        MPI_Send(&b[count - 1], 1, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD);
        MPI_Send(&d[count - 1], 1, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD);
      }
    }
  }

  float lastb = b[count - 1];
  float lastd = d[count - 1];
  for (int i = 0; i < count - 1; i++) {
    d[i] += -c[i] / lastb * lastd;
    c[i] = 0.0;
  }

  float *x = malloc(count * sizeof(float));
  for (int i = 0; i < count; i++) {
    x[i] = d[i] / b[i];
  }

  float *X = NULL;

  if (rank == 0) {
    X = malloc(n * sizeof(float));
  }

  MPI_Gatherv(x, count, MPI_FLOAT, X, sendcounts, displs, MPI_FLOAT, 0,
              MPI_COMM_WORLD);

  free(x);

  if (rank == 0) {
    free(sendcounts);
  }

  free(a);
  free(b);
  free(c);
  free(d);

  MPI_Finalize();

  if (rank == 0) {
    printvec(X, n);

    free(X);
  }

  return 0;
}
