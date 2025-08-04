clear; clc; close all;
%%

g1 = tf(0.1564, [1, 0.03101]);

Ts = 1;

g1d = c2d(g1, Ts);

[b, a] = tfdata(g1d, 'v');

c = 1;
delta = [1, -1];

a_ = conv(a,delta);

N = 70;

R = eye(N);
Q = eye(N);

A = [-a_(2:end)', [1; 0]];
B = b';
D = -a_(2:end)';

H = [1, 0];

G = zeros(N);
G_ = zeros(N,1);
F = zeros(N,length(a));
E = zeros(N,1);

for i = 1:N
    F(i,:) = H*A^(i);
    E(i) = H*A^(i-1)*D;
    G_(i) = H*A^(i-1)*B;
end

for i = 1:N
    j = i-1;
    G(i:end,i) = G_(1:end-j);
end

K = ((G'*Q*G+R)\G')*Q;
K = K(1,:);

Kr = sum(K);
KF = K*F;
KE = K*E;