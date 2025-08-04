g1 = tf(1, [0.7, 1]);   % y1/u
g2 = tf(5, [0.3, 1]);   % y2/u

Ts = 0.03;

g1d = c2d(g1, Ts);
g2d = c2d(g2, Ts);

set(g1d, 'Variable', 'z^-1');
set(g2d, 'Variable', 'z^-1');

sys = [g1d; g2d];

[num1, den1] = tfdata(g1d, 'v');
[num2, den2] = tfdata(g2d, 'v');

delta = [1 -1];
A1_ = conv(den1, delta);  % a1(q⁻¹) * (1 - q⁻¹)
A2_ = conv(den2, delta);  % a2(q⁻¹) * (1 - q⁻¹)

A1 = [ -A1_(2), 1;
       -A1_(3), 0 ];

A2 = [ -A2_(2), 1;
       -A2_(3), 0 ];

A = blkdiag(A1, A2);  % 4×4 matrix

B1 = num1';  % already discrete-time numerator
B2 = num2';

B = [
    B1(1);
    B1(2);
    B2(1);
    B2(2)
];

c = [1; 0];  % observer structure
D1 = c - A1_(2:3)';
D2 = c - A2_(2:3)';

D = [D1;
     D2];

H1 = [1, 0];  % y1 = x1
H2 = [1, 0];  % y2 = x3
H = blkdiag(H1, H2);  % 2×4 matrix

N = 10;
ny = 2;  % number of outputs
nu = 1;  % number of inputs

G = zeros(ny*N, nu*N);
F = zeros(ny*N, length(A));
E = zeros(ny*N, 1);

for i = 1:N
    F((i-1)*ny+1:i*ny, :) = H * A^i;
    E((i-1)*ny+1:i*ny) = H * A^(i-1) * D;

    for j = 1:i
        G((i-1)*ny+1:i*ny, (j-1)*nu+1:j*nu) = H * A^(i-j) * B;
    end
end

Q = eye(ny*N);      % Output tracking cost
R = 0.05 * eye(N);  % Input increment cost

K = ((G'*Q*G + R) \ G') * Q;  % (N×2N)
K = K(1,:);                   % First control move (1×2N)

Kr = K * ones(ny*N, ny);  % (1×2), gain for reference
KF = K * F;               % (1×4)
KE = K * E;               % (1×1)
