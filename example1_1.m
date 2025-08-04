clear; clc;

% plant definition
a = [1, -0.8];
b = [0, 0.4, 0.6];

% disturbance model
c = 1;
delta = [1, -1];

a_ = conv(a,delta);

N = 3;

R = 0.8*eye(N);
R(1) = 0;
Q = eye(N);

A = [-a_(2:end)', [1; 0]];
B = b(2:3)';
D = -a_(2:end)';

H = [1, 0];

G  = zeros(N);
G_ = zeros(N,1);
F  = zeros(N,length(a));
E  = zeros(N,1);

for i = 1:N
    F(i,:) = H*A^(i);
    E(i)   = H*A^(i-1)*D;
    G_(i) = H*A^(i-1)*B;
end 

for i=1:N
    j = i-1;
    G(i:end,i) = G_(1:end-j);
end

K = (G'*G+R)\G';
K = K(1,:);

Kr = sum(K);
KF = K*F;
KE = K*E;

syms z

V = KE + (KF - KE*H)*((z*eye(2) - A + D*H)\D);
P = b(2:3)/a;