function [grad,dummy] = LSS_Shur(t,a,rho,sigma,beta)
alfa = 10;
% dt = t(2)-t(1);
N = size(t,1)-1;
m = size(a,2);
B = sparse(m*N,(N+1)*m);
RHS = sparse(m*N,1);
dt = t(2:end)-t(1:end-1);
dtFrac = dt/(t(end)-t(1));
wb = 0.5*([dtFrac;0]+[0;dtFrac]);
for i=2:m
    wb = [wb,0.5*([dtFrac;0]+[0;dtFrac])];
end
wb = 1./wb;
wb = reshape(wb',1,[]); %?!

% wBinv = sparse(diag(wb));
wBinv = sparse([1:length(wb)],[1:length(wb)],wb);

uMid = 0.5*(a(2:end,:)+a(1:end-1,:));
dudt = (a(2:end,:)-a(1:end-1,:))./dt;

for i=1:N
    dfdu = [-sigma, sigma, 0; rho-uMid(i,3), -1, -uMid(i,1); uMid(i,2), uMid(i,1), -beta];
    dfds = [0, uMid(i,1), 0]';
    E = -eye(m)/dt(i) - 0.5*dfdu;
    G = eye(m)/dt(i) - 0.5*dfdu;
    b = dfds;
    
    i_start = 1+(i-1)*m;
    j_start = 1+(i-1)*m;
    i_end = i_start+(m-1);
    j_end = j_start+(m-1);
    B([i_start:i_end],[j_start:j_end]) = E;
    B([i_start:i_end],[j_end+1:j_end+1+(m-1)]) = G;
    RHS([i_start:i_end])=b;
    %     disp(i);
    %     [whySingular, ia, ic] = unique(A,'rows','stable');
end
S = B*wBinv*B';
if alfa>0
    E = blkdiag(dudt(1,:)');
    for i=2:size(dudt,1)
        E = sparse(blkdiag(E,dudt(i,:)'));
    end
    we = dtFrac*alfa*alfa;
%     wEinv = sparse(diag(1./we));
    wEinv = sparse([1:length(we)],[1:length(we)],1./we);
    S = S + E*wEinv*E';
end
condest(S)
% beep;



%%TANGENT
w = S\RHS;
v = wBinv*(B'*w);
v = reshape(v,[],size(a,1))';
if alfa>0
    eta = wEinv*(E'*w);
end
J0 = uMid(:,3); 

grad1 = mean(sum(([0,0,1].*v)')) - mean((J0-mean(J0)).*eta);
grad2 = 0;
grad = grad1 + grad2;
dummy = 0;
end