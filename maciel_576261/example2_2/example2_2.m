clear; clc;

p1 = [-0.2, 1];
p2 = [1, -1, 0];
Ps = tf(p1,p2);

Ts = 0.01;

z = tf('z',Ts);

Pz = c2d(Ps, Ts, 'zoh');

[b, a] = tfdata(Pz,'v');

% disturbance model
c1 = [1, -0.9];
c = conv(c1, c1);
c = conv(c,c1);

delta = [1, -1];

a_ = conv(a,delta);

N = 70;

R = 0*eye(N);
% R(1) = 0;
Q = eye(N);

A = [-a_(2:end)', [1; 0; 0], [0; 1; 0]];
% B = [b(2:end)'; 0];
B = b';
D = (c-a_)';
D = D(2:end);

H = [1, 0, 0];

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

G = G(:,1);
R = 0;

K = ((G'*Q*G+R)\G')*Q;
K = K(1,:);

Kr = sum(K)
KF = K*F
KE = K*E

% Kr = 6.2047
% KF = [29725, 28791, 27872]
% Kr = 6.2047;

syms z

% V = KE + (KF - KE*H)*((z*eye(2) - A + D*H)\D);
% P = b(2:3)/a;