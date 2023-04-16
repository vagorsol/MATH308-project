tempfunc = @(t,Y) kuznetov_fun(t,Y);
% make t larger until it looks like the paper's

max_time = 10^6;

% original conditions
%{
[t0, Y0] = ode45(tempfunc, [0, max_time], [10^6 10^6]); 
plot(Y0(:,1), Y0(:,2))
hold on
%}

% different inital conditions
[t2, Y2] = ode45(tempfunc, [0, max_time], [2*10^6 4.5*10^8]); 
plot(Y2(:,1), Y2(:,2))

% graph formatting things
title('ODE of Kuznetsov et al. (1994) model')
xlabel('E (cells)')
ylabel('T (cells)')

% Y1(end,:
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