clear; close all; clc;

% Time parameters
t0 = 0;          % initial value
tf = 370*22;       % final value
tspan = [t0 tf];

% Initial conditions
V0  = 0.2;
Pw0 = 0.05;
SOD0 = 1;
X0  = [V0; Pw0; SOD0];

% Resolution
[t, X] = ode15s(@three_ode, tspan, X0);
V  = X(:,1);
Pw = X(:,2);
SOD = X(:,3);

% Plot, function of time
figure(1);
plot(t/365, X(:,1), 'g', 'LineWidth', 2); hold on;
plot(t/365, X(:,2), 'b', 'LineWidth', 2);
ylim([0 1]);
xlabel('years');
legend('V','Pw')
title('Phosphorus and Vegetation as functions of time')

figure(2);
plot(Pw, V, 'g', 'LineWidth', 1.5);
xlabel('Pw (g m^{-3})');
ylabel('V (fraction de couverture)');
title('Phase diagram');
grid on;
