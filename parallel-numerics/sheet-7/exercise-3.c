#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

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

struct CoordinateForm {
  int nnz;
  int *x;
  int *y;
  double *values;
};

struct CSR {
  int nnz;
  int *rows;
  int *cols;
  double *values;
};

struct CSC {
  int nnz;
  int *rows;
  int *cols;
  double *values;
};

struct CSRdiag {
  int nnz;
  int *cols;
  double *values;
};

struct Diagonals {
  int nd;
  int *diags;
  double *values;
};

struct RectRow {
  int nl;
  int *cols;
  double *values;
};

struct JaggedDiag {
  int nnz;
  int ndiags;
  int *P;
  int *diags;
  int *cols;
  double *values;
};

double *zeros(int n) {
  double *x = malloc(n * sizeof(double));

  for (int i = 0; i < n; i++) {
    x[i] = 0.0;
  }

  return x;
}

double *MAdd(void *a, double *b, int n) {
  double *A = a;
  double *x = zeros(n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      x[i] += A[i * n + j] * b[j];
    }
  }

  return x;
}

double *CFAdd(void *a, double *b, int n) {
  struct CoordinateForm *A = a;
  double *x = zeros(n);

  for (int i = 0; i < A->nnz; i++) {
    x[A->x[i]] += A->values[i] * b[A->y[i]];
  }

  return x;
}

double *CSRAdd(void *a, double *b, int n) {
  struct CSR *A = a;
  double *x = zeros(n);

  for (int i = 0; i < n; i++) {
    for (int j = A->rows[i]; j < A->rows[i + 1]; j++) {
      x[i] += A->values[j] * b[A->cols[j]];
    }
  }

  return x;
}

double *CSCAdd(void *a, double *b, int n) {
  struct CSC *A = a;
  double *x = zeros(n);

  for (int j = 0; j < n; j++) {
    for (int i = A->cols[j]; i < A->cols[j + 1]; i++) {
      x[A->rows[i]] += A->values[i] * b[j];
    }
  }

  return x;
}

double *CSRdiagAdd(void *a, double *b, int n) {
  struct CSRdiag *A = a;
  double *x = malloc(n * sizeof(double));

  for (int i = 0; i < n; i++) {
    x[i] = A->values[i] * b[i];

    for (int j = A->cols[i]; j < A->cols[i + 1]; j++) {
      x[i] += A->values[j] * b[A->cols[j]];
    }
  }

  return x;
}

double *DiagonalsAdd(void *a, double *b, int n) {
  struct Diagonals *A = a;
  double *x = zeros(n);

  for (int i = 0; i < A->nd; i++) {
    int d = A->diags[i];

    if (d <= 0) {
      for (int j = 0; j < n + d; j++) {
        x[-d + j] += A->values[i * n + j - d] * b[j];
      }
    } else {
      for (int j = d; j < n; j++) {
        x[j - d] += A->values[i * n + j - d] * b[j];
      }
    }
  }

  return x;
}

double *RectRowAdd(void *a, double *b, int n) {
  struct RectRow *A = a;
  double *x = zeros(n);

  for (int i = 0; i < n; i++) {
    int base = i * A->nl;

    for (int j = 0; j < A->nl; j++) {
      int ind = base + j;

      x[i] += A->values[ind] * b[A->cols[ind]];
    }
  }

  return x;
}

double *JaggedAdd(void *a, double *b, int n) {
  struct JaggedDiag *A = a;
  double *x = zeros(n);

  for (int i = 0; i < A->ndiags; i++) {
    for (int j = A->diags[i]; j < A->diags[i + 1]; j++) {
      int row = j - A->diags[i];

      x[A->P[row]] += A->values[j] * b[A->cols[j]];
    }
  }

  return x;
}

const long NANO = 10 * 10 * 10 * 10 * 10 * 10 * 10 * 10 * 10;

double nanotime() {
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);

  return (double)ts.tv_sec + (double)ts.tv_nsec / (double)NANO;
}

void bm(char *msg, double *(*f)(void *, double *, int), void *A, double *b,
        int n) {
  int N = 10;
  double time = 0.0;

  for (int i = 0; i < N; i++) {
    double start = nanotime();

    free(f(A, b, n));

    time += nanotime() - start;
  }

  printf(msg, time / N);
}

int main(int argc, char **argv) {
  int n = 4000;
  double p = 0.1;
  double *A = malloc(n * n * sizeof(double));
  double *b = malloc(n * sizeof(double));
  struct CoordinateForm ACF;
  struct CSR ACSR;
  struct CSC ACSC;
  struct CSRdiag ACSRdiag;
  struct Diagonals ADiagonals;
  struct RectRow ARectRow;
  struct JaggedDiag AJaggedDiag;

  srand(5);

  for (int i = 0; i < n; i++) {
    b[i] = (float)rand() / RAND_MAX;
  }

  int nnz = 0;
  for (int i = 0; i < n * n; i++) {
    if (rand() <= p * RAND_MAX) {
      nnz += 1;
      A[i] = (float)rand() / RAND_MAX;
    }
  }

  bm("Full: %fs\n", MAdd, A, b, n);

  ACF.nnz = nnz;
  ACF.x = malloc(nnz * sizeof(int));
  ACF.y = malloc(nnz * sizeof(int));
  ACF.values = malloc(nnz * sizeof(double));

  int k = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double a = A[i * n + j];

      if (a != 0.0) {
        ACF.x[k] = i;
        ACF.y[k] = j;
        ACF.values[k] = a;
        k++;
      }
    }
  }

  bm("Coordinate Form: %fs\n", CFAdd, &ACF, b, n);

  free(ACF.x);
  free(ACF.y);
  free(ACF.values);

  ACSR.nnz = nnz;
  ACSR.rows = malloc((n + 1) * sizeof(int));
  ACSR.cols = malloc(nnz * sizeof(int));
  ACSR.values = malloc(nnz * sizeof(double));

  k = 0;
  for (int i = 0; i < n; i++) {
    ACSR.rows[i] = k;

    for (int j = 0; j < n; j++) {
      double a = A[i * n + j];

      if (a != 0.0) {
        ACSR.cols[k] = j;
        ACSR.values[k] = a;

        k++;
      }
    }
  }

  ACSR.rows[n] = nnz;

  bm("Compressed Sparse Row: %fs\n", CSRAdd, &ACSR, b, n);

  free(ACSR.rows);
  free(ACSR.cols);
  free(ACSR.values);

  ACSC.nnz = nnz;
  ACSC.rows = malloc(nnz * sizeof(int));
  ACSC.cols = malloc((n + 1) * sizeof(int));
  ACSC.values = malloc(nnz * sizeof(double));

  k = 0;
  for (int j = 0; j < n; j++) {
    ACSC.cols[j] = k;

    for (int i = 0; i < n; i++) {
      double a = A[i * n + j];

      if (a != 0.0) {
        ACSC.rows[k] = i;
        ACSC.values[k] = a;

        k++;
      }
    }
  }

  ACSC.cols[n] = nnz;

  bm("Compressed Sparse Column: %fs\n", CSCAdd, &ACSC, b, n);

  free(ACSC.rows);
  free(ACSC.cols);
  free(ACSC.values);

  int nnzdiag = 0;
  for (int i = 0; i < n; i++) {
    if (A[i * n + i] != 0.0) {
      nnzdiag += 1;
    }
  }

  ACSRdiag.nnz = nnz;
  ACSRdiag.cols = malloc((nnz + n - nnzdiag + 1) * sizeof(int));
  ACSRdiag.values = malloc((nnz + n - nnzdiag + 1) * sizeof(double));

  for (int i = 0; i < n; i++) {
    ACSRdiag.values[i] = A[i * n + i];
  }

  ACSRdiag.cols[n] = nnz + n - nnzdiag + 1;

  k = n + 1;
  for (int i = 0; i < n; i++) {
    ACSRdiag.cols[i] = k;

    for (int j = 0; j < n; j++) {
      if (i == j) {
        continue;
      }

      double a = A[i * n + j];

      if (a != 0.0) {
        ACSRdiag.values[k] = a;
        ACSRdiag.cols[k] = j;

        k += 1;
      }
    }
  }

  bm("CSR with Extracted Diagonal: %fs\n", CSRdiagAdd, &ACSRdiag, b, n);

  free(ACSRdiag.cols);
  free(ACSRdiag.values);

  int ndiags = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n - i; j++) {
      if (A[(i + j) * n + j] != 0.0) {
        ndiags += 1;
        break;
      }
    }
  }

  for (int j = 1; j < n; j++) {
    for (int i = 0; i < n - j; i++) {
      if (A[i * n + i + j] != 0.0) {
        ndiags += 1;
        break;
      }
    }
  }

  ADiagonals.nd = ndiags;
  ADiagonals.diags = malloc(ndiags * sizeof(int));
  ADiagonals.values = malloc(ndiags * n * sizeof(double));

  k = 0;

  for (int i = n - 1; i >= 0; i--) {
    int isdiag = 0;

    for (int j = 0; j < n - i; j++) {
      double a = A[(i + j) * n + j];

      if (a != 0.0) {
        isdiag = 1;
        ADiagonals.values[k * n + i + j] = a;
      }
    }

    if (isdiag) {
      ADiagonals.diags[k] = -i;
      k += 1;
    }
  }

  for (int j = 1; j < n; j++) {
    int isdiag = 0;

    for (int i = 0; i < n - j; i++) {
      double a = A[i * n + i + j];

      if (a != 0.0) {
        isdiag = 1;
        ADiagonals.values[k * n + i] = a;
      }
    }

    if (isdiag) {
      ADiagonals.diags[k] = j;
      k += 1;
    }
  }

  bm("Diagonal-Wise Storage: %fs\n", DiagonalsAdd, &ADiagonals, b, n);

  free(ADiagonals.diags);
  free(ADiagonals.values);

  int maxnnz = 0;

  for (int i = 0; i < n; i++) {
    int innz = 0;

    for (int j = 0; j < n; j++) {
      if (A[i * n + j] != 0.0) {
        innz += 1;
      }
    }

    if (innz > maxnnz) {
      maxnnz = innz;
    }
  }

  ARectRow.nl = maxnnz;
  ARectRow.cols = malloc(n * maxnnz * sizeof(int));
  ARectRow.values = malloc(n * maxnnz * sizeof(double));

  for (int i = 0; i < n; i++) {
    int k = 0;
    int base = i * ARectRow.nl;

    for (int j = 0; j < n; j++) {
      double a = A[i * n + j];

      if (a != 0.0) {
        ARectRow.cols[base + k] = j;
        ARectRow.values[base + k] = a;
        k += 1;
      }
    }

    for (int j = k; j < ARectRow.nl; j++) {
      ARectRow.cols[base + j] = 0;
      ARectRow.values[base + j] = 0.0;
    }
  }

  bm("Rectangular Row Storage: %fs\n", RectRowAdd, &ARectRow, b, n);

  free(ARectRow.cols);
  free(ARectRow.values);

  AJaggedDiag.nnz = nnz;
  AJaggedDiag.ndiags = maxnnz;
  AJaggedDiag.P = malloc(n * sizeof(int));
  AJaggedDiag.diags = malloc((maxnnz + 1) * sizeof(int));
  AJaggedDiag.cols = malloc(nnz * sizeof(int));
  AJaggedDiag.values = malloc(nnz * sizeof(double));

  for (int i = 0; i < n; i++) {
    AJaggedDiag.P[i] = i;
  }

  int *numnzs = malloc(n * sizeof(int));

  for (int i = 0; i < n; i++) {
    int nzs = 0;

    for (int j = 0; j < n; j++) {
      if (A[i * n + j] != 0.0) {
        nzs += 1;
      }
    }

    numnzs[i] = nzs;
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n - i - 1; j++) {
      if (numnzs[AJaggedDiag.P[j]] < numnzs[AJaggedDiag.P[j + 1]]) {
        int tmp = AJaggedDiag.P[j + 1];
        AJaggedDiag.P[j + 1] = AJaggedDiag.P[j];
        AJaggedDiag.P[j] = tmp;
      }
    }
  }

  free(numnzs);

  int *columns = malloc(n * sizeof(int));

  for (int i = 0; i < n; i++) {
    columns[i] = 0;
  }

  k = 0;
  int diaglength = n;
  for (int d = 0; d < maxnnz; d++) {
    AJaggedDiag.diags[d] = k;

    for (int i = 0; i < diaglength; i++) {
      int row = AJaggedDiag.P[i];

      while (columns[row] < n && A[row * n + columns[row]] == 0.0) {
        columns[row] += 1;
      }

      if (columns[row] >= n) {
        diaglength = i;
        break;
      } else {
        AJaggedDiag.cols[k] = columns[row];
        AJaggedDiag.values[k] = A[row * n + columns[row]];
        columns[row] += 1;
        k += 1;
      }
    }
  }

  AJaggedDiag.diags[maxnnz] = nnz;

  free(columns);

  bm("Jagged Diagonal: %fs\n", JaggedAdd, &AJaggedDiag, b, n);

  free(AJaggedDiag.P);
  free(AJaggedDiag.diags);
  free(AJaggedDiag.cols);
  free(AJaggedDiag.values);

  free(b);
  free(A);

  return 0;
}
