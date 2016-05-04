f = imread('fish_saltpepper.png');
m = size(f, 1);
n = size(f, 2);
f = im2double(f(:));
u = f;
lambda = 1;
epsilon = 0.1;

% Construct finite difference gradient operator
em = ones(m, 1);
en = ones(n, 1);
Dx = spdiags([-1 * en, en], [0, -1], n, n);
Dy = spdiags([-1 * em, em], [0, 1], m, m);
DX = kron(Dx', speye(m));
DY = kron(speye(n), Dy);
D = cat(1, kron(eye(3), DX), kron(eye(3), DY));

terminate = false;
atol = 150;
N = 0;

% Loop until gradient is sufficiently small
while ~terminate
    N = N + 1;
    
    % Some kind of decreasing step size
    tau = 1 / (8 * N);
    
    % Subgradient descent step
    Du = D * u;
    T = sqrt(sum(reshape(Du, n * m, 6).^2, 2));
    T = repmat(T, 6, 1);
    nonzero = T ~= 0;
    zeroind = T == 0;
    T(nonzero) = Du(nonzero) ./ T(nonzero);
    T(zeroind) = rand(nnz(zeroind), 1);
    grad = lambda / 2 * sign(u - f) + D' * T;
    u = u - tau * grad;
    
    normg = abs(norm(grad));
    fprintf('Norm: %f\n', normg);
    terminate = normg < atol;
end

fprintf('Number of iterations: %d\n', N);

u = reshape(u, m, n, 3);
imshow(u);
imwrite(u, 'denoised.png');