clc; clear; close all;

% state-space matrices of the model
A = [-0.011698, 0, 0.011698
      0, -0.02681, 0.011698
      0.011698, 0.011698, -0.023396];

B = [64.9351, 0
     0, 64.9351
     0, 0];

C = eye(3);

sys = ss(A,B,C,0);
G = tf(sys);

% transfer functions
G11 = G(1,1);   % Q1 to h1
G12 = G(1,2);   % Q2 to h1
G21 = G(2,1);   % Q1 to h2
G22 = G(2,2);   % Q2 to h2

G2 = [G11, G12
      G21, G22];

n = length(G2);
%% system identification

np = [1, 1
      1, 1];

nz = [0, 0
      0, 0];

name = 'g';

for i=1:2
    for j=1:2
        [y, t] = step(G2(i,j));
        u = ones(size(t));
        Ts = t(2) - t(1);
        data = iddata(y,u,Ts);
        p = np(i,j); z = nz(i,j);
        sys_tf = tfest(data, p, z);
        figure; compare(data, sys_tf);
        
        idxi = num2str(i);
        idxj = num2str(j);
        
        var = [name, [idxi,idxj]];
        assignin('base', var, sys_tf);
    end
end

g11 = tf(42.21, [1, 0.002745]);
g12 = tf(10.30, [1, 0.002378]);
g21 = tf(10.30, [1, 0.002378]);
g22 = tf(19.84, [1, 0.004738]);