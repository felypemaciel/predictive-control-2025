close all;

RespOpt = RespConfig('Amplitude',0.3e-4);

[y11, t11] = step(G11, RespOpt);
[y12, t12] = step(G12, RespOpt);
[y21, t21] = step(G21, RespOpt);
[y22, t22] = step(G22, RespOpt);

[y11s, t11s] = step(g11, RespOpt);
[y12s, t12s] = step(g12, RespOpt);
[y21s, t21s] = step(g21, RespOpt);
[y22s, t22s] = step(g22, RespOpt);

lw = 1.5;
ls = 18;
ts = 14;

fig1 = figure; 
plot(t11, y11, 'LineWidth',lw); hold on;
plot(t11s, y11s, 'LineWidth', lw, 'LineStyle','--', 'Color','k');
ylabel('Level (m)', 'FontSize', ls);
xlabel('Time (sec)', 'FontSize', ls);
xlim([0 2000]);
ax = gca;
ax.FontSize = ts; 
legend('Original', 'Simplified', 'Location','best');

fig2 = figure; 
plot(t12, y12, 'LineWidth',lw); hold on;
plot(t12s, y12s, 'LineWidth', lw, 'LineStyle','--', 'Color','k');
ylabel('Level (m)', 'FontSize', ls);
xlabel('Time (sec)', 'FontSize', ls);
xlim([0 2500]);
ax = gca;
ax.FontSize = ts; 
legend('Original', 'Simplified', 'Location','best');

fig3 = figure; 
plot(t21, y21, 'LineWidth',lw); hold on;
plot(t12s, y12s, 'LineWidth', lw, 'LineStyle','--', 'Color','k');
ylabel('Level (m)', 'FontSize', ls);
xlabel('Time (sec)', 'FontSize', ls);
xlim([0 2500]);
ax = gca;
ax.FontSize = ts; 
legend('Original', 'Simplified', 'Location','best');

fig4 = figure; 
plot(t22, y22, 'LineWidth',lw); hold on;
plot(t12s, y12s, 'LineWidth', lw, 'LineStyle','--', 'Color','k');
ylabel('Level (m)', 'FontSize', ls);
xlabel('Time (sec)', 'FontSize', ls);
xlim([0 1600]);
ax = gca;
ax.FontSize = ts; 
legend('Original', 'Simplified', 'Location','best');

path = 'C:\Users\room2\OneDrive\Documentos\GitHub\predictive-control-2025\paper\figures';

figu = [fig1, fig2, fig3, fig4];

for i=1:4
    idx = num2str(i);
    exportgraphics(figu(i),strcat(path,'\step',idx,'.pdf'),'Resolution',300);
end
