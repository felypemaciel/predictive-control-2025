clc; clear; close all;

% state-space matrices of the model
A = readmatrix('A2.csv');
B = readmatrix('B2.csv');
C = readmatrix('C.csv');

sys = ss(A,B,C,0);
G = tf(sys);

% transfer functions
G11 = G(1,1);   % Qin to VA
G12 = G(1,2);   % Qw to VA
G21 = G(7,1);   % Qin to COD
G22 = G(7,2);   % Qw to COD
G31 = G(6,1);   % Qin to Qgas
G32 = G(6,2);   % Qw to Qgas


G2 = [G11, G12
      G21, G22];

n = length(G2);
%% system identification

np = [1, 1
      3, 3];

nz = [0, 0
      2, 1];

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
