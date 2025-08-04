clear; clc;

A   = [1, -0.8];    % matrix A
nA  = length(A);    % A dimension
D   = [1 -1];       % Delta
A_  = conv(D, A);   % matrix Ã

N = 3;              % horizon

B = [0.4 0.6];      % matrix B
L = 0.8;            % lambda
R = L*eye(N);       % matrix R

p1 = [1, zeros(1,nA)];  % 1st polynomial
p2 = A_;                % 2nd polynomial

F   = [];   % matrix F
Gp  = [];

for n = 1:N
    [q,r]   = polydiv(p1,p2);       % polynomial divison
    F       = [F; r(2:end)];        % matrix F
    p1      = [r(2:end),0];         % update 1st polynomial
    E(n)    = q;                    % update matrix E
    aux     = conv(E,B);            % auxiliary matrix
    Gp      = [Gp; aux(n+1:end)];   % matrix Gp
end

Ga = aux(1:n)';  % matrix Ga
G  = zeros(N,N); % matrix G

% obtain triangular version of Ga
for i=1:N
    j = i-1;
    G(i:end,i) = Ga(1:end-j);
end

aux = inv(G'*G + R)*G';
k   = aux(1,:);

kr  = sum(k);
kGp = [0, k*Gp];
kF  = k*F;

Ts = 0.1;

delta = tf(1,[1 -1],Ts,'Variable','z^-1');