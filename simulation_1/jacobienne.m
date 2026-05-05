%% Paramètres de l'article
values = num2cell([0.2 30 1 0.2 4 0.05 0 0.07 475 3]);
[C_p, gamma_e, h_e, h_V, p_e, ~, r_b, r_v, tau, z] = deal(values{:});

%%Plage de P_w,in
P_w_in_min = 0.01; 
P_w_in_max = 0.1;   
n_points = 200;
P_w_in_list = linspace(P_w_in_min, P_w_in_max, n_points);

trace_list = zeros(1, n_points);
det_list = zeros(1, n_points);
V_eq_list = zeros(1, n_points);
P_eq_list = zeros(1, n_points);

for i = 1:n_points
    P_w_in = P_w_in_list(i);
    
    %calcul équilibre
    V_eq = 1 / (1 + (gamma_e * P_w_in / h_e)^p_e);
    P_eq = (z * P_w_in / C_p) * (1 + V_eq/h_V);
    
    %Variables
    P_w = P_eq/z * C_p * h_V/(h_V + V_eq);
    E = gamma_e * P_w;
    
    %Dérivées (plus simple pour calculer)
    dE_dV = -gamma_e * P_eq/z * C_p * h_V/(h_V + V_eq)^2;
    dE_dP = gamma_e * C_p/z * h_V/(h_V + V_eq);
    
    %Jacobienne
    A = r_v * (1 - 2*V_eq*(h_e^p_e + E^p_e)/h_e^p_e) ...
        - r_v * V_eq^2 * p_e*E^(p_e-1)/h_e^p_e * dE_dV;
    
    B = -r_v * V_eq^2 * p_e*E^(p_e-1)/h_e^p_e * dE_dP;
    
    C = 1/tau * C_p * h_V/(h_V + V_eq)^2 * P_eq;
    
    D = -1/tau * C_p * h_V/(h_V + V_eq);
    
    %Tr et Det
    trace_list(i) = A + D;
    det_list(i) = A*D - B*C;
end

figure(1);
hold on;
plot(P_w_in_list, det_list, 'r-', 'LineWidth', 2, 'DisplayName', 'det(J)');
xlabel('P_{w,in} (g m^{-3})', 'FontSize', 11);
ylabel('Det(J)', 'FontSize', 11);
title('Déterminant de P_{w,in}', 'FontSize', 12);
grid on;
hold off;

figure(2);
hold on
plot(P_w_in_list, trace_list, 'b-', 'LineWidth', 2, 'DisplayName', 'tr(J)');
xlabel('P_{w,in} (g m^{-3})', 'FontSize', 11);
ylabel('Tr(J)', 'FontSize', 11);
title('Trace de P_{w,in}', 'FontSize', 12);
grid on;
hold off;
