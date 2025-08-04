    clear; clc

% continuous transfer functions
g11 = tf(1, [0.7, 1]);
g12 = tf(5, [0.3, 1]);
g21 = tf(1, [0.5, 1]);
g22 = tf(2, [0.4, 1]);

Ts = 0.03;

g11d = c2d(g11, Ts);
g12d = c2d(g12, Ts);
g21d = c2d(g21, Ts);
g22d = c2d(g22, Ts);

set(g11d, 'Variable', 'z^-1');
set(g12d, 'Variable', 'z^-1');
set(g21d, 'Variable', 'z^-1');
set(g22d, 'Variable', 'z^-1');

[num11, den11] = tfdata(g11d, 'v');
[num12, den12] = tfdata(g12d, 'v');
[num21, den21] = tfdata(g21d, 'v');
[num22, den22] = tfdata(g22d, 'v');

A1_ = conv(conv(den11, den12), [1, -1]);
A2_ = conv(conv(den21, den22), [1, -1]);

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

c = [0; 0; 0];

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

% G(:,end-1:end) = []; % forgetting the two last columns

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

kr = zeros(2);

kr(1,1) = sum(K(1,1:2:end));
kr(1,2) = sum(K(1,2:2:end));
kr(2,1) = sum(K(2,1:2:end));
kr(2,2) = sum(K(2,2:2:end));