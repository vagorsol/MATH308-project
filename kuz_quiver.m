% variables
s = 13000;
d = 0.0412;
p = 0.1245;
a = 0.18;
g = 2.019.*10^7;
m = 3.422.*(10^(-10));
b = 2.*10^(-9);
n = 1.101.*10^(-7);

tempfunc = @(t,Y) kuznetov_fun(t,Y);
max_time = 10^6;

% add effector cells
[t0, Y0] = ode45(tempfunc, [0, max_time], [2*10^6 4.5*10^8]); 
plot(Y0(:,1), Y0(:,2))
hold on 

x = 0:10^2:10^6;
y = 0:10^2:10^6;

u = zeros(2, length(10^4));
v = zeros(2, length(10^4));

for i = 1:length(x)
    for j = 1:length(y)
        u(i,j) = s - d.*x(i) + p.*x(i).*(y(j)./(g+y(j))) - m.*x(i).*y(j);
        v(i,j) = a.*y(j).*(1-b.*y(j)) - n.*x(i).*y(j);;
        
        % scale arrows to same length
        c = sqrt(u(i,j)^2 + v(i,j)^2); 
        u(i,j) = .5*u(i,j)/c;
        v(i,j) = .5*v(i,j)/c;
    end
end
u = u';
v = v';
quiver(x, y, u, v)

% graph formatting things
title('ODE of Kuznetsov et al. (1994) model')
axis([0 3.5*10^6 0 5*10^8])
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