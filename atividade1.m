clear all;
close all;
clc;


a = 0.9;
b = 0.2;
N = 5;

F = zeros(N,1);
G = b*eye(N);-

for i = 1:N
    F(i) = a^i;
    if i >= 2
        for j = 1:i-1
           G(i,j) = a^(i-1)*b;
        end
    end
end

