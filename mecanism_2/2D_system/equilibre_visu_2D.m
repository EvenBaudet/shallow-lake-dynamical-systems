
c = constants();
V = linspace(0.01, 0.99, 500);

% 1. Nullcline de V (dV/dt = 0, V ≠ 0)

P_w_V = (c.hv + V) ./ (c.hv * c.gammae) .* (c.he^c.pe .* (1./V - 1)).^(1/c.pe);

% 2. Nullcline de P_w (dP_w/dt = 0)
term = (c.kv/c.l * V).^c.psod;
P_w_P = c.pwin + c.s * c.g .* term ./ (term + c.hsod^c.psod);

% Graph of nullclines
figure('Position', [100 100 800 600]);
plot(V, P_w_V, 'g-', 'LineWidth', 2); hold on;
plot(V, P_w_P, 'r-', 'LineWidth', 2);

% Search for intersections 
diff = abs(P_w_V - P_w_P);
[~, idx] = min(diff);
V_eq_est = V(idx);
P_w_eq_est = P_w_V(idx);

% Printing the intersection
if diff(idx) < 0.1  
    plot(V_eq_est, P_w_eq_est, 'bo', 'MarkerSize', 12, 'LineWidth', 2);
    text(V_eq_est, P_w_eq_est, sprintf('  (V=%.3f, P_w=%.3f)', ...
        V_eq_est, P_w_eq_est), 'FontSize', 10, 'Color', 'b');
end

xlabel('Couverture végétale V');
ylabel('Phosphore dans l''eau P_w (g/m³)');
title('Nullclines du modèle 2, 2D');
legend('dV/dt = 0 (V nullcline)', 'dP_w/dt = 0 (Pw nullcline)', ...
       'Intersection', 'Location', 'best');
grid on;
xlim([0 1]);
ylim([0 max([P_w_V, P_w_P])*1.1]);

