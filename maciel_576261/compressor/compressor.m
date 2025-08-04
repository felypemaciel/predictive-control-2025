clear; % clc

delta = [1 -1];

s = tf('s');

% continuous transfer functions
g11 = tf(0.1133, [1.783, 4.48, 1], 'IODelay', 0.715);
g12 = tf(0.9222, [2.071, 1]);
g21 = tf(0.3378, [0.361, 1.09, 1], 'IODelay', 0.299);
g22 = tf(-0.321, [0.104, 2.463, 1], 'IODelay', 0.94);

N1 = 20;
N2 = 23;
N3 = 3;
lambda = 0.8;

Max = max(N2,N3);
Min = min(N2,N3);

Ts = 0.05;

g11d = c2d(g11, Ts);
g12d = c2d(g12, Ts);
g21d = c2d(g21, Ts);
g22d = c2d(g22, Ts);

delta = tf(1, delta, Ts, 'variable', 'z^-1');
Delta = [delta, 0; 0, delta];

set(g11d, 'Variable', 'z^-1');
set(g12d, 'Variable', 'z^-1');
set(g21d, 'Variable', 'z^-1');
set(g22d, 'Variable', 'z^-1');

sys = [g11d, g12d
       g21d, g22d];

[num11, den11] = tfdata(g11d, 'v');
[num12, den12] = tfdata(g12d, 'v');
[num21, den21] = tfdata(g21d, 'v');
[num22, den22] = tfdata(g22d, 'v');

A1_ = conv(conv(den11, den12), [1, -1]);
A2_ = conv(conv(den21, den22), [1, -1]);

c = 1;

A1 = [-A1_(2), 1, 0;
      -A1_(3), 0, 1;
      -A1_(4), 0, 0];

A2 = [-A2_(2), 1, 0;
      -A2_(3), 0, 1;
      -A2_(4), 0, 0];

A = [A1, zeros(3);
    zeros(3), A2];

B11 = conv(num11, den12)';
B12 = conv(num12, den11)';
B21 = conv(num21, den22)';
B22 = conv(num22, den21)';

B = [B11(2:end), B12(2:end);
     0, 0;
     B21(2:end), B22(2:end);
     0, 0];

% syms c1 c2 c3

c = [1; 0; 0];

D1 = c - A1_(2:end)';
D2 = c - A2_(2:end)';

D = [D1, zeros(3,1);
     zeros(3,1), D2];

H1 = [1, 0, 0];
H2 = [1, 0, 0];

H = [H1, [0, 0, 0];
     [0, 0, 0], H2];

% digits(4)
% sympref('FloatingPointOutput',true);

%% building G
N = 3;
G = zeros(2*N);

for i = 1:N
    for j = 1: N
        if i >= j
            G(2*i-1:2*i,2*j-1:2*j) = H*A^(i-j)*B;
        end
    end
end

G(:,end-1:end) = 0; % forgetting the two last columns

%% building F
F = zeros(2*N);
for i=1:N
    F(2*i-1:2*i,:) = H*A^i;
end

%% building E
E = [];

for i=1:N
    E = [E; H*A^(i-1)*D];
end

l1 = 0.05; l2 = 0.05;
d1 = 1; d2 = 1;

L = repmat([l1,l2], 1, N);
d = repmat([d1,d2], 1, N);

R = diag(L);
Q = diag(d);

K_ = ((G'*Q*G + R)\G')*Q;
K = K_(1:2,:);

KF = K*F;
KE = K*E;

kr = zeros(2);

kr(1,1) = sum(K(1,1:2:end));
kr(1,2) = sum(K(1,2:2:end));
kr(2,1) = sum(K(2,1:2:end));
kr(2,2) = sum(K(2,2:2:end));