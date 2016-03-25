#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

int pow2(int b) {
  int pow = 1;

  while (b-- > 0) {
    pow *= 2;
  }

  return pow;
}

void swap(float **a, float **b) {
  float *tmp = *a;
  *a = *b;
  *b = tmp;
}

float *solve(int n, float *a, float *b, float *c, float *d) {
  int N = 1 + (int)log2((float)n);
  int bytes = n * sizeof(float);
  float *a1 = malloc(bytes);
  float *b1 = malloc(bytes);
  float *c1 = malloc(bytes);
  float *d1 = malloc(bytes);
  float *a2 = malloc(bytes);
  float *b2 = malloc(bytes);
  float *c2 = malloc(bytes);
  float *d2 = malloc(bytes);

  memcpy(a1, a, bytes);
  memcpy(b1, b, bytes);
  memcpy(c1, c, bytes);
  memcpy(d1, d, bytes);

  for (int k = 1; k <= N; k++) {
    int p2 = pow2(k - 1);

    for (int i = 0; i < n; i++) {
      float alpha = -a1[i];
      float gamma = -c1[i];
      int imp2 = i - p2;
      int ipp2 = i + p2;

      if (imp2 >= 0) {
        alpha /= b1[imp2];
      }

      if (ipp2 < n) {
        gamma /= b1[ipp2];
      }

      if (imp2 < 0 && ipp2 >= n) {
        a2[i] = 0.0;
        b2[i] = b1[i];
        c2[i] = 0.0;
        d2[i] = d1[i];
      } else if (imp2 < 0) {
        a2[i] = 0.0;
        b2[i] = b1[i] + gamma * a1[ipp2];
        c2[i] = gamma * c1[ipp2];
        d2[i] = d1[i] + gamma * d1[ipp2];
      } else if (ipp2 >= n) {
        a2[i] = alpha * a1[imp2];
        b2[i] = alpha * c1[imp2] + b1[i];
        c2[i] = 0.0;
        d2[i] = alpha * d1[imp2] + d1[i];
      } else {
        a2[i] = alpha * a1[imp2];
        b2[i] = alpha * c1[imp2] + b1[i] + gamma * a1[ipp2];
        c2[i] = gamma * c1[ipp2];
        d2[i] = alpha * d1[imp2] + d1[i] + gamma * d1[ipp2];
      }
    }

    swap(&a1, &a2);
    swap(&b1, &b2);
    swap(&c1, &c2);
    swap(&d1, &d2);
  }

  float *x = malloc(bytes);

  for (int i = 0; i < n; i++) {
    x[i] = d1[i] / b1[i];
  }

  free(a1);
  free(b1);
  free(c1);
  free(d1);
  free(a2);
  free(b2);
  free(c2);
  free(d2);

  return x;
}

int main(int argc, char **argv) {
  int n = 5;

  float a[] = {0.0, -1.0, -1.0, -1.0, -1.0};
  float b[] = {2.0, 2.0, 2.0, 2.0, 2.0};
  float c[] = {-1.0, -1.0, -1.0, -1.0, 0.0};
  float d[] = {2.5, 1.5, 2.0, 0.1, 5.0};

  float *x = solve(n, a, b, c, d);

  printvec(x, n);

  free(x);

  return 0;
}
