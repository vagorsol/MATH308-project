% variables
s = 13000;
d = 0.0412;
p = 0.1245;
a = 0.18;
g = 2.019.*10^7;
m = 3.422.*(10^(-10));
b = 2.*10^(-9);
n = 1.101.*10^(-7);

x = 0:10^2.5:10^6;
y = 0:10^2.5:10^6;

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