function predatorNpreyS(x,y)
num=1000;
dt=0.12;
u = zeros(num,1);
v = zeros(num,1);
syms m n
u(1)=x;
v(1)=y;

for k = 2:num
eqns = [m-u(k-1) == dt*u(k-1)*(n-2), n-v(k-1) == dt*n*(1-u(k-1))];


[x1,x2] = solve(eqns, [m n]);
u(k)=x1(1);
v(k)=x2(1);
end

plot(u,v,'.');
axis equal;
xlim([0,10]);
ylim([0,10]);