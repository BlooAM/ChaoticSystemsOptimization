clear
clc

sigma = 10;
beta = 8/3;
z_mean = [];
for rho=25:0.1:100
    z_mean = [z_mean,lorenz_int(sigma,rho,beta)];
end
rho = [25:0.1:100];
figure(1)
plot(rho,z_mean)
grid on
xlabel('rho')
ylabel('z_mean')

r0 = 28;
r1 = 28+(1E-60j);
z_mean0 = lorenz_int(10,r0,8/3);
z_mean1 = lorenz_int(10,r1,8/3);
d_z_mean_d_r = (z_mean1-z_mean0)/(r1-r0);


function z_mean = lorenz_int(s,r,b)
T = 100;
dt = 0.01;
x = 1;
y = 2;
z = 3;
z_int = 0;
for t=0:dt:T
    dxdt = s*(y-x);
    dydt = x*(r-z)-y;
    dzdt = x*y-b*z;
    x = x+dt*dxdt;
    y = y+dt*dydt;
    z = z+dt*dzdt;
    z_int = z_int+dt*z;
end
z_mean = z_int/T;
end