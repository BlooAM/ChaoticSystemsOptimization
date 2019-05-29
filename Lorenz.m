clear
clc

sigma = 10;
b = 8/3;
for r = 0:10
dy = @(t,y)[sigma*(y(2)-y(1));-y(1)*y(3)+r*y(1)-y(2);y(1)*y(2)-b*y(3)];
[t,trajectory] = ode45(dy,[0 100],[0 0.5 1]);
figure(1)
p = plot3(trajectory(:,1),trajectory(:,2),trajectory(:,3))
hold on
grid on
end