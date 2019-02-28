A = [1,5,7;2,10,3;3,3,-19];
b = [1,6,-1]';

x0 = [1,1,1]';
eps = 1e-5;
iter = 1;
x = x0;
while iter<10000%max(abs(A*x-b)>eps)
    x(1) = b(1) - A(1,2:end)*x(2:end);
    for i = 2:length(x)-1
       x(i) = (b(i)-A(i,1:i-1)*x(1:i-1)-A(i,i+1:end)*x(i+1:end))/A(i,i); 
    end
    x(end) = b(end) - A(end,1:end-1)*x(1:end-1);
    disp(A*x-b);
    iter = iter+1;
end