M = [
  1  1  0  0  0  0  0;
  0  0  1  1  0  0  0;
  -1 0 -1  0  1  0  0;
  0  0  0  0 -1  1  1;
  0 -1  0  0  0 -1  0;
  0  0  0 -1  0  0 -1];
z = [2 2 2 2 1 2 2]';
d = [1 0 0 0 0 -1]';

C = -M' * d / norm(d)^2;
A = M - d * d' * M / norm(d)^2;
b = zeros(size(d));
lb = zeros(size(z));
ub = z;

glpk(C, A, b, lb, ub)
