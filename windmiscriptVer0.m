%% Matlab script version of the windmi model %%
%% Ver0. Goals 
% 1.Solve connected ODEs of the windmi-RC model.
% 2. Tunable/time-varying parameters
%% %%%%%%%%%%%%%%
clc;close all;
tspan = [0 20];
y0 = 1;
[t1,y1] = ode23s(@(t,y) -10*t, tspan, y0);
plot(t1,y1,'r--');
hold on; grid on; axis tight
[t2,y2] = ode15s(@(t,y) -10*t, tspan, y0);
plot(t2,y2,'k--')
