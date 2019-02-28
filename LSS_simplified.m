function [dJds_mod,dJds] = LSS_simplified(t,a,rho,sigma,beta)
alfa = 10;
dt = t(2)-t(1);
m = size(t,1)-2;
n = size(a,2);
B = sparse(n*m,m*(n+1)+n);
RHS = sparse(n*m,1);
for i=2:m+1
    dfdu_minus = [-sigma, sigma, 0; rho-(a(i-1,3)+a(i,3))/2, -1, -(a(i-1,1)+a(i,1))/2; (a(i-1,2)+a(i,2))/2, (a(i-1,1)+a(i,1))/2, -beta];
    dfdu_plus = [-sigma, sigma, 0; rho-(a(i,3)+a(i+1,3))/2, -1, -(a(i,1)+a(i+1,1))/2; (a(i,2)+a(i+1,2))/2, (a(i,1)+a(i+1,1))/2, -beta];
    dfds_minus = [0, (a(i-1,1)+a(i,1))/2, 0]';
    dfds_plus = [0, (a(i,1)+a(i+1,1))/2, 0]';
    u_plus = ((a(i+1,:)+a(i,:))')/2;
    u_minus = ((a(i,:)+a(i-1,:))')/2;
    E = -eye(n)/dt - 0.5*dfdu_minus;
    G = eye(n)/dt - 0.5*dfdu_plus;
    f = (u_plus-u_minus)/dt;
    b = (dfds_minus+dfds_plus)/2;
    
    i_start = 1+(i-2)*n;
    j_start = 1+(i-2)*(n+1);
    i_end = i_start+(n-1);
    j_end = j_start+(n-1);
    B([i_start:i_end],[j_start:j_end]) = E;
    B([i_start:i_end],j_end+1) = f/alfa;
    B([i_start:i_end],[j_end+2:j_end+2+(n-1)]) = G;
    RHS([i_start:i_end])=b;
%     disp(i);
%     [whySingular, ia, ic] = unique(A,'rows','stable');
end
A = B*B';
x = A\RHS;
w = x;
condest(A)
beep;

ve = sparse(n*m,1);
eta = sparse(m-1,1);
dfdu_minus = [-sigma, sigma, 0; rho-(a(1,3)+a(2,3))/2, -1, -(a(1,1)+a(2,1))/2; (a(1,2)+a(2,2))/2, (a(1,1)+a(2,1))/2, -beta];
E = -eye(n)/dt - 0.5*dfdu_minus;
ve(1:3) = -E'*w(1:3);
for i=2:m
    i_start = 1+(i-2)*n;
    j_start = 1+(i-2)*(n+1);
    i_end = i_start+(n-1);
    j_end = j_start+(n-1);
    
    ve(3+i_start:3+i_end) = -B([i_start:i_end],[j_end+2:j_end+2+(n-1)])'*w(i_start:i_end)-B([i_end+1:i_end+n],[j_end+2:j_end+2+(n-1)])'*w(i_end+1:i_end+n);
    eta(i-1) = -(B([i_start:i_end],j_end+1)*alfa)'*w(i_start:i_end)/(alfa*alfa);
%     B([i_start:i_end],[j_start:j_end]) = E;
%     B([i_start:i_end],j_end+1) = f/alfa;
%     B([i_start:i_end],[j_end+2:j_end+2+(n-1)]) = G;

%     disp(i)
end

J = 0;
Jbar = 0;
for i=1:m
    Jbar = Jbar+(a(i,3)+a(i+1,3))/2;
end
Jbar = Jbar/m;
for i=2:m
    J = J+eta(i-1)*(((a(i-1,3)+a(i,3))/2+(a(i,3)+a(i+1,3))/2)/2-Jbar);
end
J = J/(m-1);
temp = 0;
for i=1:m
    temp = temp+[0,0,1]*ve(1+(i-1)*n:1+(i-1)*n+n-1);
end
temp = temp/m;
dJds = temp+J;
dJds_mod = J-temp;
disp(J-temp)
end