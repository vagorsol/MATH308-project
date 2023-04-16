% variables
alpha = 1.636;
beta = 2*10^(-3);
sigma = 0.1181;
delta = 0.3643;
mu = 0.00311;
rho = 1.131;
eta = 20.19; 

% E-nullcline function
efunc = @(y) sigma./(delta + mu.*y - ((rho.*y)./(eta + y)));

% E-nullcline func plot pt.1
y = 0:0.1:10;
plot(efunc(y), y, 'Color',[0.4940 0.1840 0.5560])
hold on

% T-nullcline (first one)
y=500:-0.1:0;
tfunc = @(y) alpha.*(1 - beta.*y);
plot(tfunc(y), y, 'Color', [0.4660 0.6740 0.1880])

% E-nullcline func plot pt.2
y=500:-0.1:220;
plot(efunc(y), y, 'Color',[0.4940 0.1840 0.5560])

% T-nullcline (second one)
x = 0:1:5;
plot(x,zeros(size(x)), 'Color', [0.4660 0.6740 0.1880]) 

% plot intersection points
labels = {'A - unstable', 'B - stable', 'C - unstable', 'D - stable'};
e = [0.3155 1.6093 0.7172 0.1825];
t = [0 8.158 280.8 442.2];
plot(e, t, 'o') 

text(e, t, labels, 'VerticalAlignment','bottom','HorizontalAlignment','left');
title('Nullclines of Kuznetsov et al. (1994) model')
legend('E-nullclines','T-nullclines')
axis([0 5 0 500])
xlabel('Effector Cells (x)')
ylabel('Tumor Cells (y)')