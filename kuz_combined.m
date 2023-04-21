tempfunc = @(t,Y) kuznetov_fun(t,Y);
max_time = 10^6;

% nullcline variables
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
y=10^6:-10^2:0;
tfunc = @(y) alpha.*(1 - beta.*y);
plot(tfunc(y), y, 'Color', [0.4660 0.6740 0.1880])

% E-nullcline func plot pt.2
y=10^6:-10^2:220;
plot(efunc(y), y, 'Color',[0.4940 0.1840 0.5560])

% T-nullcline (second one)
x = 0:10^2:10^6;
plot(x,zeros(size(x)), 'Color', [0.4660 0.6740 0.1880]) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ode system graph
[t0, Y0] = ode45(tempfunc, [0, max_time], [10^6 10^6]); 
plot(Y0(:,1), Y0(:,2))

% graph formatting things
title('ODE of Kuznetsov et al. (1994) model')
legend('E-nullclines','T-nullclines', 'Small Tumor Equilibrium', 'Large Tumor State', 'Added Effector Cells')
axis([0 10^6 0 10^6])
xlabel('E (cells)')
ylabel('T (cells)')

function v = kuznetov_fun(t,Y)
    s = 13000;
    d = 0.0412;
    p = 0.1245;
    a = 0.18;
    g = 2.019.*10^7;
    m = 3.422.*(10^(-10));
    b = 2.*10^(-9);
    n = 1.101.*10^(-7);

    E = Y(1);
    T = Y(2);
    
    v(1) = s - d.*E + p.*E.*(T./(g+T)) - m.*E.*T;
    v(2) = a.*T.*(1-b.*T) - n.*E.*T;
    
    v = v';
end
