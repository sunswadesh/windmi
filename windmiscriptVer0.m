%% Matlab script version of the windmi model %%
%% Ver0. Goals 
% 1.Solve connected ODEs of the windmi-RC model.
% 2. Tunable/time-varying parameters
%% %%%%%%%%%%%%%%
clc;close all;

% t0 = 0;
% tfinal = 15;
% % p0 = [50; 50];
% p0 = -50
% [t,p] = ode45(@lotkaODE,[t0 tfinal],p0);
% plot(t,p)
% title('Predator/Prey Populations Over Time')
% xlabel('t')
% ylabel('Population')
% legend('Prey','Predators')

ft = linspace(0,155,50);
f = ones(1,length(ft)); f(end/2:end) = 0; % Create step input

gt = linspace(0,155,50);
% g = 3*sin(gt-0.25);
g = ones(1,length(gt))*10; g(end/1.8:end) = 50; % Two step ring current recovery
y0 = [10;10];

tspan = [1 155];
ic = [10;10];
opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
[t,y] = ode45(@(t,y) myode(t,y,ft,f,gt,g), tspan, ic, opts);

plot(t,-y,'*-'); grid on; axis tight
%The ODEs.
% function dpdt = lotkaODE(t,p)        
% % LOTKA Lotka-Volterra predator-prey model
% delta = 0.02;
% beta = 0.01;
% % 
% % dpdt = [p(1) .* (1 - beta*p(2));
% %         p(2) .* (-1 + delta*p(1))];
% 
% dpdt = [-p(1) ./10];
% end
function dydt = myode(t,y,ft,f,gt,g)
f = interp1(ft,f,t); % Interpolate the data set (ft,f) at time t
g = interp1(gt,g,t); % Interpolate the data set (gt,g) at time t
dydt = [(-y(1)./g + f); y(1).*y(2)];% Evaluate ODE at time t
end