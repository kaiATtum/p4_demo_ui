
res1 = load('rho_1_e_min_2.mat');
res2 = load('rho_1_e_min_3.mat');
res3 = load('rho_5_e_min_2.mat');
res4 = load('rho_5_e_min_3.mat');

figure()
figSize=[2, 2, 6, 2];
fntsize=10;
set(0,'DefaultAxesFontSize',fntsize);
set(gcf,'Units','inches','Position',figSize);


plot(res1.results_file,'DisplayName','\rho_0 = 1e-2'); hold on;
plot(res2.results_file,'DisplayName','\rho_0 = 1e-3'); hold on;
plot(res3.results_file,'DisplayName','\rho_0 = 5e-2'); hold on;
plot(res4.results_file,'DisplayName','\rho_0 = 5e-3'); hold off;

set(gca, 'YScale', 'log')
xlabel('iterations');
ylabel('convergence error');
filename='residual';
grid on;
legend('Location','North','Orientation','horizontal');
legend show;

export_fig(gcf,filename,'-nocrop','-pdf','-r300','-painters');