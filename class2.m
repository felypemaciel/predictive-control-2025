% calculate the prediction of the disturbance model using the Diophantine
% equation for
% a = 0.9
% b = 0.2
% N = 5

% calculate F, G, and the disturbance model's predictions.

a = 0.9;
b = 0.2;
N = 5;

u = zeros(N,1);

F = zeros(N,1);
G = zeros(N,N);
U = zeros(N,1);

Gcol = zeros(N,1);

for i=1:N
    F(i) = a^i;
    Gcol(i) = a^(i-1)*b;
    U(i) = u(i);
end

for i=1:N
    j = i-1;
    G(i:end,i) = Gcol(1:end-j);
end