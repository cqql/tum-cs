s = 2;
N = 4;

% Original vector
x = zeros(N, 1);
ind = randperm(N);
x(ind(1:s)) = normrnd(3, 1, s, 1);
% x = [0; 1; 0; 1];

% DFT matrix
DFT = zeros(2*s, N);
for i = 1:2*s
    for j = 1:N
        DFT(i, j) = exp(2 * pi * 1i * (i - 1) * (j - 1) / N);
    end
end

% Measurement
y = DFT * x;
% Measurement error
% e = normrnd(0, 0.1, 2*s, 1);
e = zeros(2 * s, 1);

% Approximate p(1)..p(s)
A = zeros(s);
for i = 1:s
    A(i, :) = fliplr(y(i:(i + s - 1))');
end
q = A \ -y((s + 1):2*s);
q = q / N;

% Evaluate polynomial q on 1:N
vals = zeros(1, N);
for i = 1:N
    vals(i) = 1 + exp(2 * pi * 1i * i * (1:s) / N) * q;
end

% S is the set of indices with the smallest values
[V, I] = sort(abs(vals));
S = I(1:s);

% Invert S
S = 1 + N - S;

% Solve for x
X = zeros(N, 1);
X(S) = DFT(:, S) \ y;