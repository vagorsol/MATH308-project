tempfunc = @(t,Y) kuznetov_fun(t,Y);
max_time = 10^6;

% add effector cells
[t0, Y0] = ode45(tempfunc, [0, max_time], [2*10^6 4.5*10^8]); 
plot(Y0(:,1), Y0(:,2))
hold on 

% trajectory going to large tumor state
[t1, Y1] = ode45(tempfunc, [0, max_time], [0 0.15*10^8]); 
plot(Y1(:,1), Y1(:,2))

% dotted line (added E cells
plot([0.0017*10^8 2*10^6], [4.4728*10^8 4.5*10^8], 'LineStyle', '--')
plot(0.0017*10^8, 4.4728*10^8, 'o','Color','black')

% graph formatting things
title('ODE of Kuznetsov et al. (1994) model')
axis([0 3.5*10^6 0 5*10^8])
legend('Small Tumor Equilibrium', 'Large Tumor State', 'Added Effector Cells')
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