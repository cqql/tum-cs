f = imread('noisy_input.png');
m = size(f, 1);
n = size(f, 2);
f = im2double(f(:));
u = f;
lambda = 5;
epsilon = 0.1;

% Half of the upper bound
tau = 1 / (2 * lambda + normest(D)^2 / epsilon);

% Construct finite difference gradient operator
em = ones(m, 1);
en = ones(n, 1);
Dx = spdiags([-1 * en, en], [0, -1], n, n);
Dy = spdiags([-1 * em, em], [0, 1], m, m);
DX = kron(Dx', speye(m));
DY = kron(speye(n), Dy);
D = cat(1, kron(eye(3), DX), kron(eye(3), DY));

terminate = false;
atol = 3;
N = 0;

% Loop until gradient is sufficiently small
while ~terminate
    % Gradient descent step
    Du = D * u;
    dh = Du ./ sqrt(Du.^2 + epsilon^2);
    grad = lambda * (u - f) + D' * dh;
    u = u - tau * grad;
    
    normg = abs(norm(grad));
    fprintf('Norm: %f\n', normg);
    terminate = normg < atol;
    N = N + 1;
end

fprintf('Number of iterations: %d\n', N);

u = reshape(u, m, n, 3);
imshow(u);
imwrite(u, 'denoised.png');