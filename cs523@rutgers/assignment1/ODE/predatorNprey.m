function predatorNprey(x,y)
n=100;
dt=0.12;
u = zeros(n,1);
v = zeros(n,1);

u(1)=x;
v(1)=y;

for k = 2:n
   u(k)=dt*u(k-1)*(v(k-1)-2)+u(k-1);
   v(k)=dt*v(k-1)*(1-u(k-1))+v(k-1);
end

plot(u,v,'.');
axis equal;
xlim([0,10]);
ylim([0,10]);