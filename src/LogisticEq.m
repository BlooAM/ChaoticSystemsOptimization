%0<r<4, 0<x0<1
r = 3.9;
x0 = 0.9;
nIter = 1000;

x = [x0];
for i=2:nIter
    x = [x,logistic(r,x(end))];
end


figure(1)
plot(x)
grid on
drawnow

function x_n = logistic(r,x)
x_n = r*x*(1-x);
end