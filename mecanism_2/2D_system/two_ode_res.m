clear; close all; clc;

% --- Time parameters
t0 = 0;          % initial value
tf = 16*365;       % final value
tspan = [t0 tf];

% --- Initial conditions
V0  = 0.2;
Pw0 = 0.05;
X0  = [V0; Pw0];

% --- Resolution
[t, X] = ode15s(@two_ode, tspan, X0);
V  = X(:,1);
Pw = X(:,2);

% --- Plot
figure;
plot(t/365, X(:,1), 'g', 'LineWidth', 2); hold on;
plot(t/365, X(:,2), 'b', 'LineWidth', 2);
ylim([0 1]);    
xlabel('Years');
title('Phosphorus and Vegetation as functions of time, 2D-system')
legend('V','Pw')