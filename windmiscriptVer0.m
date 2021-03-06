%% Matlab script version of the windmi model %%
%% Ver0. Goals 
% 1.Solve connected ODEs of the windmi-RC model.
% 2. Tunable/time-varying parameters
%% %%%%%%%%%%%%%%
clc;clear;%close all;

etm = 1*24*3600;
vswt = linspace(0,etm,etm);
vsw = zeros(1,length(vswt)); vsw(1:25000) = 1e3; % Create step input

rct = linspace(0,etm,etm);
tau_rc = ones(1,length(rct))*2*3600; %tau_rc(end/2:end) = 50*3600; % Two step ring current recovery

tspan = [1 etm];
ic = [0; 0; 0; 0; 0; 0; 0; 0];
% ic = [0; 0; 0];
opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
[t,y] = ode45(@(t,y) myode(t,y,vswt,vsw,tau_rc), tspan, ic, opts);

vrstr = {'I'; 'V'; 'p'; 'Kp'; 'I1'; 'Vi'; 'I2'; 'Wrc'};
for i = [1,8]
    figure(i)
    subplot(211);plot(vswt/3600,vsw); grid on; axis tight; hold on
    subplot(212);plot(t/3600,-y(:,i),'*-'); grid on; axis tight; hold on
    title(vrstr{i})
end

figure(2);
funt = (1 + tanh((y(:,1) - 350e4)./125))./2;
plot(t/3600,funt);grid on;hold on; axis tight
% plot(t,)

function dydt = myode(t,y,vswt,vsw,tau_rc)
%% Windmi constants
    R_Earth = 6380e3;  % radius if earth in meters
    L = 100; del_I = 1.25e+05;
    M = 2/3;    C = 48000;    Sigma = 7.8;    Omega = 10000*R_Earth^3;
    mu0 = 4.2e-09;    Ic = 1.78e+8;    alpha = 7.550246e+11;
    taup = 10*60;    tauE = 0.5*3600;    L1 = 19;    C1 = 800;
    SigmaI = 3;    beta = 1;    alpha_F = 100;    R_prc = 0.1;
%     tau_rc = 12*3600;    
    L2 = 8;    R_A2 =0.3;    B_tr = 5e-9;
    A_eff = 2*(R_Earth^2);    L_y = 5*R_Earth;    Alph = 4.3;
%% %%%%%%%%%%%%% 
%% The 8 unknowns are: I, V, p, Kp, I1, Vi, I2, Wrc. y(1-8)          
% windmi ODEs without the mutual inductance terms
Ic = 350e4; del_I = 125;
vsw = interp1(vswt,vsw,t); % Interpolate the data set (vswt,vsw) at time t
tau_rc = interp1(vswt,tau_rc,t); % Interpolate the data set (gt,tau_rc) at time t
dydt = [(vsw - y(2)./L);
    (y(1) - y(5) - alpha.*sqrt(y(3)) - Sigma.*y(2))./C; % Ips is derived from plasma sheet pressure
    (Sigma.*y(2).^2./Omega - mu0.*y(3).*sqrt(y(4)).*(1 + tanh((y(1) - Ic)./del_I))./2 ...,
    - y(3).*y(2).*A_eff./(Omega.*B_tr.*L_y) - 1.5*y(3)./tauE)*(2/3); % No substorm
    alpha.*sqrt(y(3)).*y(2) - y(4)./taup;
    (y(2) - y(6))./L1;
    (y(5) - y(7) -SigmaI.*y(6))./C1;
    (y(6) - (R_prc + R_A2).*y(7))./L2;
    R_prc.*y(7).^2 + y(3).*y(2).*A_eff./(B_tr.*L_y) - y(8)./tau_rc];
end





