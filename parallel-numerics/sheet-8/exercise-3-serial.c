#include <stdio.h>
#include <stdlib.h>

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
  int kstop = 20;
  double omega = 1.5;
  int n = 4;
  double *A = malloc(n * n * sizeof(A[0]));
  double *b = malloc(n * sizeof(b[0]));

  srand(5);

  /* A[0] = 10.0; */
  /* A[1] = 0; */
  /* A[2] = 1; */
  /* A[3] = 0.5; */
  /* A[4] = 7; */
  /* A[5] = 1; */
  /* A[6] = 1; */
  /* A[7] = 0; */
  /* A[8] = 6; */

  /* b[0] = 21; */
  /* b[1] = 9; */
  /* b[2] = 8; */

  A[0] = 4.0;
  A[1] = -1.0;

  for (int i = 1; i < n - 1; i++) {
      A[i * n + i - 1] = -1.0;
      A[i * n + i] = 4.0;
      A[i * n + i + 1] = -1.0;
  }

  A[n * n - 2] = -1.0;
  A[n * n - 1] = 4.0;

  for (int i = 0; i < n; i++) {
    b[i] = (double)rand() / RAND_MAX;
  }

  double *alpha = malloc(n * sizeof(alpha[0]));
  double *c = malloc(n * n * sizeof(c[0]));
  double *x = malloc(n * sizeof(x[0]));

  for (int i = 0; i < n; i++) {
    alpha[i] = omega / A[i * n + i];
    x[i] = 0.0;
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i == j) {
        c[i * n + j] = (1.0 - omega) * A[i * n + i] / omega;
      } else {
        c[i * n + j] = -A[i * n + j];
      }
    }
  }

  for (int k = 0; k < kstop; k++) {
    for (int i = 0; i < n; i++) {
      double s = b[i];

      for (int j = 0; j < n; j++) {
        s += c[i * n + j] * x[j];
      }

      x[i] = alpha[i] * s;
    }
  }

  printvec(x, n);

  free(x);
  free(c);
  free(alpha);

  free(b);
  free(A);

  return 0;
}
