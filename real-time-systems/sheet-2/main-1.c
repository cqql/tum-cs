#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

double F(double x) {
  return (pow(x, 1./5.) + sin(x * M_PI / 180.0)) / (x - 1./3.) * pow(x, 8.) / (-1. * cos(x * M_PI / 180.0));
}

void print_random() {
  printf("%d\n", rand() % 1000000);
}

int main() {
  printf("%e\n", F(10));

  srand(time(NULL));
  print_random();

  return 0;
}
