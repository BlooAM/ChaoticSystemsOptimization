sigma = 10;
b = 8/3;
r = 10;
rho = [];
dy={}; y={}; y={};
nIter = 3;
for i=1:nIter
    r = 10*i;
    rho(i) = r;
    dy{i} = @(t,y)[sigma*(y(2)-y(1));
                -y(1)*y(3)+r*y(1)-y(2);
                y(1)*y(2)-b*y(3)];
    % rozwi¹zanie uk³adu        
    [t{i},y{i}] = ode45(dy{i},[0 10],[1 0 28]);
end
% rysowanie wyniku
figure(1)
for i=1:nIter
    plot3(y{i}(:,1),y{i}(:,2),y{i}(:,3),'DisplayName',['rho = ',num2str(rho(i))])
    grid on
    legend()
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold on
    
end