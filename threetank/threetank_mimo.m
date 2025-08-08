clear all; clc; % close all;

S = 0.0154;     % tank cross sectional area (m2)
Sp = 5E-5;      % inter tank cross sectional area (m2)

mu = 0.5;       % outflow coefficients
mu20 = 0.675;

qmax = 1.2E-4;  % maximum flowrate (m3/s)
lmax = 0.6;     % maximum level (m)

g = 9.8;       % gravity (m/s2)
%% Model
% Process
G1 = tf(42.21, [1, 0.002745]);
G2 = tf(10.30, [1, 0.002378]);
G3 = tf(10.30, [1, 0.002378]);
G4 = tf(19.84, [1, 0.004738]);
GG = [G1 G2; G3 G4];

% Sampling time
Ts = 0.5;

GGd = c2d(GG,Ts,'zoh');  % Discretization

%% Controller tuning parameters
% Prediction horizon
N = 100;
% Control horizon
Nu = 2;

% Weight values
R = eye(2*Nu);
Q = 0.5*eye(2*N);

% Constraints
du_min = -0.5e-4;   % m3/s
du_max = 0.5e-4;    % m3/s
u_min = 0;          % m3/s
u_max = 1.2e-4;     % m3/s
y_min = 0;          % m
y_max = 0.6;        % m

%% CARIMA model
delta = [1 -1];
Delta = tf(1,delta,Ts,'Variable','z^-1'); 

[b1,a1] = tfdata(GGd(1,1),'v');
[b2,a2] = tfdata(GGd(1,2),'v');
[b3,a3] = tfdata(GGd(2,1),'v');
[b4,a4] = tfdata(GGd(2,2),'v');

c1  = conv([1 0],[1 0]);
c2  = conv([1 0],[1 0]);

a1t = conv(delta,conv(a1,a2));
a2t = conv(delta,conv(a3,a4));
b1t = conv(b1,a2);
b2t = conv(b2,a1);
b3t = conv(b3,a4);
b4t = conv(b4,a3);

% adjusting the vectors to have the same length
Max = max([length(a1t)-1,length(a2t)-1,length(b1t)-1,length(b2t)-1,length(b3t)-1,length(b4t)-1,length(c1)-1,length(c2)-1]);
a1t = [a1t zeros(1,Max-(length(a1t)-1))];
a2t = [a2t zeros(1,Max-(length(a2t)-1))];
b1t = [b1t zeros(1,Max-(length(b1t)-1))];
b2t = [b2t zeros(1,Max-(length(b2t)-1))];
b3t = [b3t zeros(1,Max-(length(b3t)-1))];
b4t = [b4t zeros(1,Max-(length(b4t)-1))];
c1 = [c1 zeros(1,Max-(length(c1)-1))];
c2 = [c2 zeros(1,Max-(length(c2)-1))];

%% State-Space formulation observable cannonical form
A1 = [-a1t(2:end)' [eye(length(a1t)-2); zeros(1,length(a1t)-2)]];
A2 = [-a2t(2:end)' [eye(length(a2t)-2); zeros(1,length(a2t)-2)]];
A = [A1 zeros(size(A1)); zeros(size(A2)) A2];
B1=[b1t(2:end)]';
B2=[b2t(2:end)]';
B3=[b3t(2:end)]';
B4=[b4t(2:end)]';
B=[B1 B2; B3 B4];
D1 = c1(2:end)'-a1t(2:end)';
D2 = c2(2:end)'-a2t(2:end)';
D = [D1 zeros(size(D1)); zeros(size(D2)) D2];
H = [1 zeros(1,Max-1)];
H = [H zeros(size(H)); zeros(size(H)) H];

%% Prediction Matrix
G = zeros(2*N,2*Nu);
for i=1:N
    for j=1:Nu
        if i>=j
            G(2*i-1:2*i,2*j-1:2*j) = H*A^(i-j)*B;
        end
    end
end

F = zeros(N * size(H,1), size(A,2));
for i = 1:N
    row_idx = (i-1)*size(H,1)+1 : i*size(H,1);
    F(row_idx, :) = H * A^i;
end

E = zeros(N * size(H,1), size(D,2));
for i = 1:N
    row_idx = (i-1)*size(H,1)+1 : i*size(H,1);
    E(row_idx, :) = H * A^(i-1) * D;
end

%% Simulation parameters

% Time parameters
Tfinal = 600;   % seconds
n_samples = ceil(Tfinal/Ts);

% References
r1 = 0.2;
r2 = 0.15;

%% GPC
% Quadratic equation HH
HH = 2*(G'*Q*G + R);

Acon = [ 1  zeros(1,2*Nu-1)
        -1  zeros(1,2*Nu-1)
         1  0 zeros(1,2*Nu-2) 
         0  1 zeros(1,2*Nu-2)
        -1  0 zeros(1,2*Nu-2) 
         0 -1 zeros(1,2*Nu-2) 
         G
        -G];

% Inicial conditions
f = zeros(2*N,1);
r = zeros(2*N,1);
y0 = zeros(2,1);
du0 = zeros(2,1);
u0 = zeros(2,1);
x0 = zeros(length(A), 1);
e0 = zeros(2,1);
[~,z01] = filter(b1(2:end),a1,0);
[~,z02] = filter(b2(2:end),a2,0);
[~,z03] = filter(b3(2:end),a3,0);
[~,z04] = filter(b4(2:end),a4,0);

du = zeros(2,n_samples);

tic;
% return
for k=1:n_samples
    y(:,k) = y0;
    x(:,k) = A*x0 + B*du0 + D*e0;
    e(:,k) = y(:,k) - H*x(:,k);
    
    % free response
    for i=1:2*N
       f(i,1) = F(i,:)*x(:,k) + E(i,:)*e(:,k); 
    end

    % references
    if k > 150/Ts
        r1 = 0.3;
    end
    if k > 250/Ts
        r1 = 0.4;
    end
    if k > 350/Ts
        r1 = 0.35;
    end

    % r2 = r1/2;
    r = repmat([r1; r2], ceil(2*N/2), 1);
    r = r(1:2*N);
    ref(:,k) = r;
    
    % quadratic equation bb
    bb = 2*(f-r)'*Q*G;

    Bcon = [du_max
           -du_min
           (u_max-u0)
          -(u_min-u0)
            y_max-f
           -y_min+f];

    % QP solver
    opt = optimoptions('quadprog','Display','off');
    % Acon = []; Bcon = []; 
    sol = quadprog(HH,bb,Acon,Bcon,[],[],[],[],[],opt);

    du(:,k) = sol(1:2);
    % disp(du(:,k));
    u(:,k) = u0 + du(:,k);
    disp(u(:,k));

    [y1,z01] = filter(b1(2:end),a1,u(1,k),z01);
    [y2,z02] = filter(b2(2:end),a2,u(2,k),z02);
    Y1 = y1+y2;

    [y3,z03] = filter(b3(2:end),a3,u(1,k),z03);
    [y4,z04] = filter(b4(2:end),a4,u(2,k),z04);
    Y2 = y3+y4;

    du0 = du(:,k);
    u0 = u(:,k);
    x0 = x(:,k);
    e0 = e(:,k);
    y0 = [Y1; Y2];
end

toc;

% return
%% Figure
t = 0:Ts:Tfinal-Ts;

figure
subplot(2,2,1)
hold on
plot(t,y(1,:), 'linewidth', 1)
plot(t,ref(1,:), '--k', 'linewidth', 1)
xlabel('t (sec)')
ylabel('h_1 (m)')

subplot(2,2,2);
hold on
plot(t,y(2,:), 'linewidth', 1)
plot(t,ref(2,:), '--k', 'linewidth', 1)
xlabel('t (sec)')
ylabel('h_2 (m)')

subplot(2,2,3)
plot(t,u(1,:), 'linewidth', 1)
xlabel('t (sec)')
ylabel('Q_1 (m^3/s)')

subplot(2,2,4)
plot(t,u(2,:), 'linewidth', 1)
xlabel('t (sec)')
ylabel('Q_2 (m^3/s)')

%% export data
data = [ref(1:2,:)', y', u'];
matrix_name = ['data',num2str(N),'.mat'];
save(matrix_name,'data');