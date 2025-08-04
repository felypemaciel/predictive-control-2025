A = [1.8 1; -0.8 0];
B = [0.4; 0.6];
D = [1.8; -0.8];
C = 0;
H = [1 0];

E = [];
F = [];
G = H*B*eye(N);

N = 3;

for i = 1:N
    % E(i) = H*A^(i-1)*D;
    E(i,1) = H*A^(i-1)*D;
    F(i,:) = H*A^i;
    if i >= 2
        for j = 1:i-1
            G(i,j) = H*A^(j)*B;
        end
    end
end

R = 0.8*eye(3);

k = inv(G'*G + R)*G';
k = k(1,:);

kr = sum(k);

kF = k*F;

kE = k*E;

z=tf('z',Ts);
b = [0.4 0.6];
a = [1 -0.8];
sys= tf(b,a,Ts,'Variable','q^-1','IODelay',1);