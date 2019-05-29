function x = LSS_optimization(alfa)
nSec = 131;
dt = 0.005;
sigma = 10;
beta = 8/3;

%Init condition
x0 = 0.1;
y0 = 0.1; %2
z0 = 0.1; %3

%Design parameter
rhoStart = 90;
rhoEnd = 90;
interval = 10;

%Containers for results
AV = [];

for rho=rhoStart:interval:rhoEnd
    f = @(t,a) [-sigma*a(1) + sigma*a(2); rho*a(1) - a(2) - a(1)*a(3); -beta*a(3) + a(1)*a(2)];
    [t,a] = ode45(f,[0:dt:nSec],[x0 y0 z0]);
    aav = [0,0,0];
    
    for i=1:size(a,1)-1
        dt = t(i+1)-t(i);
        aav =  aav + dt*(a(i,:)+a(i+1,:))/2 ;
    end
    aav = aav/(t(end)-t(1));
    AV = [AV; aav];
end
% alfa = 10;
dt = t(2)-t(1);
m = size(t,1)-1;
n = size(a,2);
B = sparse(n*(m-1),(m-1)*(n+1)+n);
RHS = sparse(n*(m-1),1);
for i=2:size(t,1)-1
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
% A = A(1:end-3,1:end-3);
% RHS = RHS(1:end-3);
x = condest(A);
end

% x0 = 1;
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
% bottomBound = 0;
% topBound = []
% [alfa_f,feval] = fmincon(@LSS_optimization,x0,[],[],[],[],bottomBound,topBound,[],options)