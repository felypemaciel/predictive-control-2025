clc; clear; close all;

a1 = [1, -0.99];
a2 = [1, -0.47];
a = conv(a1,a2);

delta = [1, -1];

N  = 6;
Nu = N;

a_  = conv(a, delta);
b1_ = conv([0, 0, 0, 0, 0.1], [1, -0.47]);

b2  = conv([1, -0.99], [0, -0.5]);
b2_ = conv(b2,delta);

c = 1;

A1 = [-a_(2:end), zeros(1,3)];
A2 = eye(N-1);
A = [A1', [A2; zeros(1,5)]];

B1 = b1_';
B2 = [b2_'; zeros(2,1)];

D = A1';

H = zeros(1,N);
H(1) = 1;

G1 = zeros(N, Nu);
G2 = zeros(N, Nu);

F = zeros(N);
E = zeros(N,1);

G1col = G1(:,1);
G2col = G1(:,2);

for i = 1:N
    G1col(i) = H*A^(i-1)*B1;
    G2col(i) = H*A^(i-1)*B2;

    F(i,:) = H*A^i;
    E(i) = H*A^(i-1)*D;
end

for i = 1:Nu
    j = i-1;
    G1(i:end,i) = G1col(1:end-j);
    G2(i:end,i) = G2col(1:end-j);
end