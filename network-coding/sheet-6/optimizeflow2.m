M = [
  1  1  0  0  0  0  0;
  0  0  1  1  0  0  0;
  -1 0 -1  0  1  0  0;
  0  0  0  0 -1  1  1;
  0 -1  0  0  0 -1  0;
  0  0  0 -1  0  0 -1];
N = [M zeros(size(M)); zeros(size(M)) M];
z = [2 2 2 2 1 2 2]';
d1 = [1 0 0 0 -1 0]';
d2 = [0 1 0 0 0 -1]';
D = [d1 zeros(size(d1)); zeros(size(d2)) d2];
DTDinv = [(1 / norm(d1)^2) 0; 0 (1 / norm(d2)^2)];

C = -(ones(1, 2) * DTDinv * D' * N)';
B = N - D * DTDinv * D' * N;
b = zeros(size(d1, 1) + size(d2, 1), 1);

# Encode Ax <= z in the constraint coefficients
I = eye(size(M, 2));
A = [I I];
ctype = [repmat("S", [size(B, 1), 1]);
         repmat("U", [size(A, 1), 1])];
B = [B; A];
b = [b; z];

[xopt, fmin, errnum, extra] = glpk(C, B, b, [], [], ctype)
