clear
clc
%ERGODICITY - independence of initial condition !!!
nSec = 100; %10
dt = 0.01;
perturbFlag = 1;
LSSFlag = 0;
sigma = 10;
beta = 8/3;

%Init condition
x0 = 1; %1
y0 = 0; %0
z0 = 28; %28

%Design parameter
rhoStart = 10;
rhoEnd = 100;
interval = 0.01;

%Containers for results
AV = [];
sensitivityLSS = [];
TESTAV = [];
Z = [];
sensitivity = [];

for rho=rhoStart:interval:rhoEnd
    %Objective function
    f = @(t,a) [-sigma*a(1) + sigma*a(2); rho*a(1) - a(2) - a(1)*a(3); -beta*a(3) + a(1)*a(2)];
%     [t,a] = ode45(f,[0 nSec],[x0 y0 z0]);
    [t,a] = ode45(f,[0:dt:nSec],[x0 y0 z0]);
    aav = [0,0,0];
    
    if LSSFlag==1
        %LSS
%         [dJds_mod,dJds] = LSS_simplified(t,a,rho,sigma,beta);
        [dJds_mod,dJds] = LSS_Shur(t,a,rho,sigma,beta);
        sensitivityLSS = [sensitivityLSS;dJds_mod,dJds];
    end
    
    for i=1:size(a,1)-1
        dt = t(i+1)-t(i);
        aav =  aav + dt*(a(i,:)+a(i+1,:))/2 ;
    end
    aav = aav/(t(end)-t(1));
    AV = [AV; aav];
    if perturbFlag == 1
        %Perturbed objective function
        perturb = (1E-60j);
        f = @(t,a) [-sigma*a(1) + sigma*a(2); (rho+perturb)*a(1) - a(2) - a(1)*a(3); -beta*a(3) + a(1)*a(2)];
%         [t,a] = ode45(f,[0 nSec],([x0 y0 z0]+perturb*v(1:3)'));
         [t_im,a_im] = ode45(f,[0:dt:nSec],([x0 y0 z0])); %perturb*v(1:3)'
        
        %Calculate sensitivity of objective function
        sensitivity = [sensitivity;imag(sum(a_im(:,3))*dt/(t_im(end)-t_im(1)))/abs(perturb)];
    end
    disp(rho)
%     figure(7)
%     plot3(a(:,1),a(:,2),a(:,3))
%     grid on
%     drawnow
end
% RESULT = [AV(:,3),sensitivity]
rho=[rhoStart:interval:rhoEnd];
%Plot result
figure(1)
plot(rho,AV(:,3))
grid on
title('Objective function')
xlabel('rho')


if perturbFlag==1 && LSSFlag==1
figure(2)
% plot()
semilogy(rho,sensitivity,rho,sensitivityLSS(:,1))
legend('Complex var method','LSS')
grid on
title('Objective function sensitivity')
xlabel('rho')
end
if LSSFlag==1
figure(3)
plot(rho,sensitivityLSS(:,1))
legend('LSS')
grid on
title('Objective function sensitivity')
ylim([-5,5])
xlabel('rho')
end
if perturbFlag==1
figure(4)
% plot(rho,sensitivity)
semilogy(rho,sensitivity)
legend('Complex var method')
grid on
title('Objective function sensitivity')
xlabel('rho')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [dJds_mod, dJds] = LSS_KKT(t,a,rho,sigma,beta)
% alfa = 1;
% dt = t(2)-t(1);
% m = size(t,1)-1;
% n = size(a,2);
% B = sparse(n*(m-1),(m-1)*(n+1)+n);
% RHS = sparse(n*(m-1),1);
% for i=2:m
%     dfdu_minus = [-sigma, sigma, 0; rho-(a(i-1,3)+a(i,3))/2, -1, -(a(i-1,1)+a(i,1))/2; (a(i-1,2)+a(i,2))/2, (a(i-1,1)+a(i,1))/2, -beta];
%     dfdu_plus = [-sigma, sigma, 0; rho-(a(i,3)+a(i+1,3))/2, -1, -(a(i,1)+a(i+1,1))/2; (a(i,2)+a(i+1,2))/2, (a(i,1)+a(i+1,1))/2, -beta];
%     dfds_minus = [0, (a(i-1,1)+a(i,1))/2, 0]';
%     dfds_plus = [0, (a(i,1)+a(i+1,1))/2, 0]';
%     u_plus = ((a(i+1,:)+a(i,:))')/2;
%     u_minus = ((a(i,:)+a(i-1,:))')/2;
%     E = -eye(n)/dt - 0.5*dfdu_minus;
%     G = eye(n)/dt - 0.5*dfdu_plus;
%     f = (u_plus-u_minus)/dt;
%     b = -(dfds_minus+dfds_plus)/2;
%     
%     i_start = 1+(i-2)*n;
%     j_start = 1+(i-2)*(n+1);
%     i_end = i_start+(n-1);
%     j_end = j_start+(n-1);
%     B([i_start:i_end],[j_start:j_end]) = E;
%     B([i_start:i_end],j_end+1) = f/alfa;
%     B([i_start:i_end],[j_end+2:j_end+2+(n-1)]) = G;
%     RHS([i_start:i_end])=b;
% %     disp(i);
% %     [whySingular, ia, ic] = unique(A,'rows','stable');
% end
% A = speye((n+1)*(m-1)+n,(n+1)*(m-1)+n);
% for i=4:4:size(A,1)
%     A(i,i) = alfa*alfa;
% end
% RHS = [sparse((n+1)*(m-1)+n,1);RHS];
% A = [A,B';B,sparse(n*(m-1),n*(m-1))];
% condest(A)
% beep;
% x = A\RHS;
% % w = x;
% ve = x(1:(n+1)*(m-1)+n);
% ve([4:4:(n+1)*(m-1)+n]) = [];
% eta = x(4:4:(n+1)*(m-1)+n);
% w = x((n+1)*(m-1)+n+1:end);
% 
% J = 0;
% Jbar = 0;
% for i=1:m-1
%     Jbar = Jbar+(a(i,3)+a(i+1,3))/2;
% end
% Jbar = Jbar/m;
% for i=2:m-1
%     J = J+eta(i-1)*(((a(i-1,3)+a(i,3))/2+(a(i,3)+a(i+1,3))/2)/2-Jbar);
% end
% J = J/(m-1);
% temp = 0;
% for i=1:m-1
%     temp = temp+[0,0,1]*ve(1+(i-1)*n:1+(i-1)*n+n-1);
% end
% temp = temp/m;
% dJds = temp+J;
% disp(J-temp)
% dJds = J-temp;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [dJds_mod,dJds] = LSS_simplified(t,a,rho,sigma,beta)
% alfa = 1;
% dt = t(2)-t(1);
% m = size(t,1)-2;
% n = size(a,2);
% B = sparse(n*m,m*(n+1)+n);
% RHS = sparse(n*m,1);
% for i=2:m+1
%     dfdu_minus = [-sigma, sigma, 0; rho-(a(i-1,3)+a(i,3))/2, -1, -(a(i-1,1)+a(i,1))/2; (a(i-1,2)+a(i,2))/2, (a(i-1,1)+a(i,1))/2, -beta];
%     dfdu_plus = [-sigma, sigma, 0; rho-(a(i,3)+a(i+1,3))/2, -1, -(a(i,1)+a(i+1,1))/2; (a(i,2)+a(i+1,2))/2, (a(i,1)+a(i+1,1))/2, -beta];
%     dfds_minus = [0, (a(i-1,1)+a(i,1))/2, 0]';
%     dfds_plus = [0, (a(i,1)+a(i+1,1))/2, 0]';
%     u_plus = ((a(i+1,:)+a(i,:))')/2;
%     u_minus = ((a(i,:)+a(i-1,:))')/2;
%     E = -eye(n)/dt - 0.5*dfdu_minus;
%     G = eye(n)/dt - 0.5*dfdu_plus;
%     f = (u_plus-u_minus)/dt;
%     b = (dfds_minus+dfds_plus)/2;
%     
%     i_start = 1+(i-2)*n;
%     j_start = 1+(i-2)*(n+1);
%     i_end = i_start+(n-1);
%     j_end = j_start+(n-1);
%     B([i_start:i_end],[j_start:j_end]) = E;
%     B([i_start:i_end],j_end+1) = f/alfa;
%     B([i_start:i_end],[j_end+2:j_end+2+(n-1)]) = G;
%     RHS([i_start:i_end])=b;
% %     disp(i);
% %     [whySingular, ia, ic] = unique(A,'rows','stable');
% end
% A = B*B';
% x = A\RHS;
% w = x;
% condest(A)
% beep;
% 
% ve = sparse(n*m,1);
% eta = sparse(m-1,1);
% dfdu_minus = [-sigma, sigma, 0; rho-(a(1,3)+a(2,3))/2, -1, -(a(1,1)+a(2,1))/2; (a(1,2)+a(2,2))/2, (a(1,1)+a(2,1))/2, -beta];
% E = -eye(n)/dt - 0.5*dfdu_minus;
% ve(1:3) = -E'*w(1:3);
% for i=2:m
%     i_start = 1+(i-2)*n;
%     j_start = 1+(i-2)*(n+1);
%     i_end = i_start+(n-1);
%     j_end = j_start+(n-1);
%     
%     ve(3+i_start:3+i_end) = -B([i_start:i_end],[j_end+2:j_end+2+(n-1)])'*w(i_start:i_end)-B([i_end+1:i_end+n],[j_end+2:j_end+2+(n-1)])'*w(i_end+1:i_end+n);
%     eta(i-1) = -(B([i_start:i_end],j_end+1)*alfa)'*w(i_start:i_end)/(alfa*alfa);
% %     B([i_start:i_end],[j_start:j_end]) = E;
% %     B([i_start:i_end],j_end+1) = f/alfa;
% %     B([i_start:i_end],[j_end+2:j_end+2+(n-1)]) = G;
% 
% %     disp(i)
% end
% 
% J = 0;
% Jbar = 0;
% for i=1:m
%     Jbar = Jbar+(a(i,3)+a(i+1,3))/2;
% end
% Jbar = Jbar/m;
% for i=2:m
%     J = J+eta(i-1)*(((a(i-1,3)+a(i,3))/2+(a(i,3)+a(i+1,3))/2)/2-Jbar);
% end
% J = J/(m-1);
% temp = 0;
% for i=1:m
%     temp = temp+[0,0,1]*ve(1+(i-1)*n:1+(i-1)*n+n-1);
% end
% temp = temp/m;
% dJds = temp+J;
% dJds_mod = J-temp;
% disp(J-temp)
% end


%% COMMENTS
% -Cel: wrazliwosc srednich calkowych trajektorii po zmiennych
%   projektyowych
% -minimalizujemy funkcjonal (24) aby aby dostac trajektorie zblizona do U
%   i przeksztalcenie czasu postaci tau = t + C
% -trajektoria jest zblizona, warunek poczakowy przesuniety - > z zal ergodycznosci
%   to nie szkodzi na srednie calkowe po duzych czasach,
%   dostajemy gladkie srednie czasowe w zaleznosci od parametrow projektowych
% -nachylenie tych srednich calkowy w zal od par. proj. zblizone do
%   "globalnego" nachylenia srednich calkowych trajektorii chaotycznych (w funkcji par. proj.)
%   czyli tego co szukamy
% -pytanie: dlaczego w (24) szukamy minimum tau i po U, a nie tylko to tau (czy U nie ma byc ustalone?)
%  innymi slowy: czy nie chodzi nam o to aby znalesc dla danego parametru
%  proj. (czyli ustalonego U) jego wplywu na wrazliwosc trajektorii
% -szukajc trajektorii: wybieramy parametry poczatkowe projektowe -> dla
%   dla nich liczymy trajekktorie (U(t)) -> zapewne jest chaotyczna ->
%   szukamy trajektorii cienia minimalizujac funkcjonal (24) szukajac
%   odpowiedniego przeksztalcenia czasu (czy to jest wykonalne - czy w odp malym otoczeniu jest
%   odpowiednio gesto tych dobrych trajektorii (uklad chaotyczny)?)
% -zakladajac ze mamy U' - ona odpowiada U dla danych par. proj., czy dla
%   innego zestawu par. proj tez bedzie jej odpowiadala ? (przekstalcenie czasu 'globalne' ?, tzn
%   srednia_calkowa(U'(tau(t)))) w funckji par. proj - gladka ?