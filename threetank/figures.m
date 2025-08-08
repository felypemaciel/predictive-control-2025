clear; close all;

t = 0:0.5:600-0.5;

N = [5 10 50 100];

load('data5.mat');
u1_5 = data(:,5);
u2_5 = data(:,6);

load('data10.mat');
u1_10 = data(:,5);
u2_10 = data(:,6);

load('data50.mat');
u1_50 = data(:,5);
u2_50 = data(:,6);

load('data100.mat');
u1_100 = data(:,5);
u2_100 = data(:,6);

fig1 = figure;
plot(t, u1_10,  'LineWidth', 1.5); hold on;
plot(t, u1_5,   'LineWidth', 1.5); 
legend('N = 10', 'N = 5');
xlabel('time (sec)');
ylabel('Q1 (m^3/s)');

fig2 = figure;
plot(t, u1_50,  'LineWidth', 1.5); hold on;
plot(t, u1_100,   'LineWidth', 1.5); 
legend('N = 50', 'N = 100');
xlabel('time (sec)');
ylabel('Q1 (m^3/s)');

fig3 = figure;
plot(t, u2_5,  'LineWidth', 1.5); hold on;
plot(t, u2_10,   'LineWidth', 1.5); 
legend('N = 5', 'N = 10');
xlabel('time (sec)');
ylabel('Q2 (m^3/s)');

fig4 = figure;
plot(t, u2_100,   'LineWidth', 1.5); hold on;
plot(t, u2_50,  'LineWidth', 1.5);
legend('N = 100', 'N = 50');

path = 'C:\Users\room2\OneDrive\Documentos\GitHub\predictive-control-2025\paper\figures';

exportgraphics(fig1, [path,'\constrained-input-Q1-1.pdf'],'Resolution',300);
exportgraphics(fig2, [path,'\constrained-input-Q1-2.pdf'],'Resolution',300);
exportgraphics(fig3, [path,'\constrained-input-Q2-1.pdf'],'Resolution',300);


% plot(t, u1_50,  'LineWidth', 1.5);
% plot(t, u1_100, 'LineWidth', 1);

% figure;
% for i = 1:length(N)
%     name_data = ['data',num2str(N(i)),'.mat'];
%     load(name_data);
%     plot(data(:,5),'LineWidth',1); hold on;
%     legend(num2str(N(i)))
% end

%% Unconstrained
path = 'C:\Users\room2\OneDrive\Documentos\GitHub\predictive-control-2025\paper\figures';
fig = figure; 
plot(t,y(1,:),'LineWidth',2); hold on;
plot(t,ref(1,:),'LineWidth',1,'LineStyle','--','Color','k');
plot(t,y(2,:),'LineWidth',2,'Color','#dd5400');
plot(t,ref(2,:)', 'LineWidth',1,'LineStyle','-.','Color','k');
legend('h_1','h_{1,ref}','h_2','h_{2,ref}','Location','best');
xlabel('time (sec)');
ylabel('level (m)');
% exportgraphics(fig,[path,'\unconstrained-control.pdf'],'Resolution',300);
exportgraphics(fig,[path,'\constrained-control.pdf'],'Resolution',300);
%%
fig = figure; 
plot(t,u(1,:),'LineWidth',2); hold on;
plot(t,u(2,:),'LineWidth',2);
legend('Q_1','Q_2','Location','best');
xlabel('time (sec)');
ylabel('Inflow (m^3/s)');
exportgraphics(fig,[path,'\unconstrained-inputs.pdf'],'Resolution',300);

%% time 

close all

x = [5 10 20 50 100];
y = [2.96 3.73 7.87 20.25 86.7];
f = fit(x', y', 'exp1');

fig = figure; % scatter(x,y,'filled')
plot(f, x, y);
xlabel('Prediction horizon (N)');
ylabel('Average running time (sec)');
legend('Data', 'Fitted curve', 'Location', 'best');
exportgraphics(fig, [path,'\running-time.pdf'],'Resolution',300);