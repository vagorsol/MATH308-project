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
labels_saddle = {'A - unstable','C - unstable'};
labels_steady = {'B - stable','D - stable'};
es = [0.3155 0.76];
ts = [0 272];

ed = [1.6093 0.1825];
td = [8.158 442.2];

plot(es, ts, 'o', 'Color', 'black') 
plot(ed, td, '.', 'Color', 'black') 

% trajectories
s = 13000;
d = 0.0412;
p = 0.1245;
a = 0.18;
g = 2.019.*10^7;
m = 3.422.*(10^(-10));
b = 2.*10^(-9);
n = 1.101.*10^(-7);

x = 0:1:5;
y = 0:100:500;

u = zeros(5);
v = zeros(5);

for i = 1:length(x)
    for j = 1:length(y)
        u(i,j) = s - d.*x(i) + p.*x(i).*(y(j)./(g+y(j))) - m.*x(i).*y(j);
        v(i,j) = a.*y(j).*(1-b.*y(j)) - n.*x(i).*y(j);
        
        % scale arrows to same length
        c = sqrt(u(i,j)^2 + v(i,j)^2); 
        u(i,j) = .5*u(i,j)/c;
        v(i,j) = .5*v(i,j)/c;
    end
end

u = u';
v = v';
quiver(x, y, u, v)

text(es, ts, labels_saddle, 'VerticalAlignment','bottom','HorizontalAlignment','left');
text(ed, td, labels_steady, 'VerticalAlignment','bottom','HorizontalAlignment','left');
title('Nullclines of Kuznetsov et al. (1994) model')
legend('E-nullclines','T-nullclines')
axis([0 5 0 500])
xlabel('Effector Cells (x)')
ylabel('Tumor Cells (y)')