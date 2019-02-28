function [dJds_mod,dJds] = LSS_adjoint(t,a,rho,sigma,beta)
alfa = 1;
dt = t(2)-t(1);
m = size(t,1)-1;
n = size(a,2);
B = sparse(n*(m-1),(m-1)*(n+1)+n);
RHS = sparse(n*(m-1),1);
for i=2:m
    dfdu_minus = [-sigma, sigma, 0; rho-(a(i-1,3)+a(i,3))/2, -1, -(a(i-1,1)+a(i,1))/2; (a(i-1,2)+a(i,2))/2, (a(i-1,1)+a(i,1))/2, -beta];
    dfdu_plus = [-sigma, sigma, 0; rho-(a(i,3)+a(i+1,3))/2, -1, -(a(i,1)+a(i+1,1))/2; (a(i,2)+a(i+1,2))/2, (a(i,1)+a(i+1,1))/2, -beta];
    u_plus = ((a(i+1,:)+a(i,:))')/2;
    u_minus = ((a(i,:)+a(i-1,:))')/2;
    E = -eye(n)/dt - 0.5*dfdu_minus;
    G = eye(n)/dt - 0.5*dfdu_plus;
    f = (u_plus-u_minus)/dt;
    b = [0,0,0]';
    
    i_start = 1+(i-2)*n;
    j_start = 1+(i-2)*(n+1);
    i_end = i_start+(n-1);
    j_end = j_start+(n-1);
    B([i_start:i_end],[j_start:j_end]) = E;
    B([i_start:i_end],j_end+1) = f;
    B([i_start:i_end],[j_end+2:j_end+2+(n-1)]) = G;
    RHS([i_start:i_end])=b;
end

Jbar = 0;
for i=1:m
    Jbar = Jbar+(a(i,3)+a(i+1,3))/2;
end
Jbar = Jbar/m;

rhs = speye((n+1)*(m-1)+n,1);
A = speye((n+1)*(m-1)+n,(n+1)*(m-1)+n);
for i=2:1:m
    A(4*(i-1),4*(i-1)) = alfa*alfa;
    rhs(4*(i-1),1) = (((a(i-1,3)+a(i,3))/2+(a(i,3)+a(i+1,3))/2)/2-Jbar)/(m-1);
    rhs(4*(i-1)-1,1) = 1/m;
end
rhs(end,1) = 1/m;
RHS = [rhs;RHS];
A = [A,B';B,sparse(n*(m-1),n*(m-1))];
condest(A)
beep;
x = A\RHS;
% w = x;
ve = x(1:(n+1)*(m-1)+n);
ve([4:4:(n+1)*(m-1)+n]) = [];
eta = x(4:4:(n+1)*(m-1)+n);
w = x((n+1)*(m-1)+n+1:end);

J = 0;
for i=2:m
    dfds_minus = [0, (a(i-1,1)+a(i,1))/2, 0]';
    dfds_plus = [0, (a(i,1)+a(i+1,1))/2, 0]';
    J = J + ((dfds_minus+dfds_plus)/2)'*w(1+(i-2)*n:1+(i-2)*n+n-1);
end

dJds = J
dJds_mod = 0;
end