clear; clc; close all;
%%

G11 = tf(0.1564, [1, 0.03101]);
G12 = tf(-0.4268, [1, 0.0326]);
G21 = tf([686.3, -3.318, 0.001205], [1, 0.07475, 0.003348, 3.586e-5]);
G22 = tf([103.9, 0.4829], [1, 0.1, 0.004264, 7.342e-5]);

delta = [1, -1];

Ts = 0.03;

G11d = c2d(G11, Ts);
G12d = c2d(G12, Ts);
G21d = c2d(G21, Ts);
G22d = c2d(G22, Ts);

delta = tf(1, delta, Ts, 'variable', 'z^-1');
Delta = [delta, 0; 0, delta];

set(G11d, 'variable', 'z^-1');
set(G12d, 'variable', 'z^-1');
set(G21d, 'variable', 'z^-1');
set(G22d, 'variable', 'z^-1');

sys = [G11d, G12d
       G21d, G22d];

[num11, den11] = tfdata(G11d,'v');
[num12, den12] = tfdata(G12d,'v');
[num21, den21] = tfdata(G21d,'v');
[num22, den22] = tfdata(G22d,'v');

A1_ = conv(conv(den11, den12), [1 -1]);
A2_ = conv(conv(den21, den22), [1 -1]);

c = 1;

A1 = [-A1_(2), 1, 0
      -A1_(3), 0, 1
      -A1_(4), 0, 0];

A2 = [-A2_(2), 1, 0
      -A2_(3), 0, 1
      -A2_(4), 0, 0];

A = [A1, zeros(3)
    zeros(3), A2];

B11 = conv(num11, den12)';
B12 = conv(num12, den11)';
B21 = conv(num21, den22)';
B22 = conv(num22, den21)';

B = [B11(2:end), B12(2:end)
    0, 0
    B21(2:end), B22(2:end)
    0, 0];

c = [1; 0; 0];

D1 = c - A1_(2:end)';
D2 = c - A2_(2:end)';

D = [D1, zeros(3,1)
    zeros(3,1), D2];

H1 = [1, 0, 0];
H2 = [1, 0, 0];

H = [H1, 0, 0, 0
     0, 0, 0, H2];

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

F = zeros(2*N);
for i=1:N
    F(2*i-1:2*i,:) = H*A^i;
end

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