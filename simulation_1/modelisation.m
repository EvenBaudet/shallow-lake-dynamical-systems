%%
values = num2cell([0.2 30 1 0.2 4 0.05 0 0.07 475 3]);
[C_p, gamma_e, h_e, h_V, p_e, P_w_in, r_b, r_v, tau, z] = deal(values{:});

T = 60*365; % en jours
dt = 0.01*365; % en jours
N = round(T/dt);

% listes pour l'intégration avec des conditions initiales proches
% du point d'équilibre

V_eq = zeros(20,N); % je rajoute le point d'équilibre après
P_eq = zeros(20,N); % pareil
P_w_eq = zeros(21,N);
E_eq = zeros(21,N);
t_eq = linspace(0,T,N);

% Point d'équilibre
V_eq_unique = 1/(1 + (gamma_e*P_w_in/h_e)^p_e);   
P_eq_unique = z*P_w_in*(h_V + V_eq_unique)/(C_p*h_V);

% calcul des courbe de "déviation" rouge et vertes
% l'évolution du lac varie beaucoup en fonction d'où on se place
% autour de cette courbe sur le diagramme P,V (une des composantes
% du gradient s'inverse)

V_deviation = linspace(0.006,0.99,200);
P_deviation_V = zeros(1,200);
P_deviation_P = zeros(1,200);

for i = 1:200
    P_deviation_V(i) = z*h_e*(h_V + V_deviation(i))*(((1 - V_deviation(i))/V_deviation(i)).^(1/p_e))/(C_p*h_V*gamma_e);
    P_deviation_P(i) = z*P_w_in*(h_V + V_deviation(i))/(C_p*h_V);
end

% Conditions initiales autour de l'équilibre

for i = 1:20
    k = find(abs(V_deviation - V_eq_unique) < 0.0025);
    V_eq(i,1) = V_deviation(k - 10 + i);
    P_eq(i,1) = P_deviation_V(k - 10 + i);
end
% rajout du point d'équilibre
V_eq = [V_eq(1:10, :); [V_eq_unique, zeros(1,N-1)]; V_eq(11:20, :)];
P_eq = [P_eq(1:10, :); [P_eq_unique, zeros(1,N-1)]; P_eq(11:20, :)];

P_w_eq(:,1) = P_eq(:,1).*C_p.*h_V./(z.*(h_V + V_eq(:,1)));
E_eq(:,1) = gamma_e*P_w_eq(:,1);

% Intégration temporelle
for i = 1:21
    for k = 2:N
        V_eq(i,k) = V_eq(i,k-1) + (1 - V_eq(i,k-1)*(h_e^p_e + E_eq(i,k-1)^p_e)/(h_e^p_e))*r_v*V_eq(i,k-1)*dt;
        P_eq(i,k) = P_eq(i,k-1) + ((P_w_in - P_w_eq(i,k-1))*z/tau - r_b*(P_eq(i,k-1) - z*P_w_eq(i,k-1))) * dt;
        P_w_eq(i,k) = P_eq(i,k)*C_p*h_V/(z*(h_V + V_eq(i,k)));
        E_eq(i,k) = gamma_e * P_w_eq(i,k);
    end
end

% Figures
figure(1)
hold on
for i = 1:21
    if i == 11
    plot(t_eq/365, squeeze(V_eq(i,:)), 'o-','LineWidth', 1, 'MarkerSize', 0.5, 'Color', [0.9, 0.1, 0.3]);
    else
    rgb = [0.3 - 0.01*i, 0.6 - 0.02*i, 1 - 0.02*i];
    plot(t_eq/365, squeeze(V_eq(i,:)), 'o-','LineWidth', 1, 'MarkerSize', 2, 'Color', rgb);
    end

end
plot(t_eq/365, squeeze(V_eq(i,:)), 'o-','LineWidth', 1, 'MarkerSize', 0.5, 'Color', rgb)
xlabel('Temps (années)');
ylabel('Végétation');
title('Végétation en fonction du temps');
grid on
hold off
    
    
% Figure Végétation vs Phosphore
figure(3)  
hold on
for i = 1:21
    rgb = [0, 0.5 - 0.02*i, 1 - 0.02*i];
    plot(squeeze(P_eq(i,:)), squeeze(V_eq(i,:)), 'o-','LineWidth', 1, 'MarkerSize', 2, 'Color', rgb);
end
plot(P_deviation_V(1:end-5), V_deviation(1:end-5), 'o-','LineWidth', 2, 'MarkerSize', 2, 'Color', [0.9, 0.1, 0.3]);
plot(P_deviation_P(28:60), V_deviation(28:60), 'o-','LineWidth', 2, 'MarkerSize', 2, 'Color', [0.3, 0.8, 0.5])
xlabel('Phosphore total (P) g.m-2'); %27:60
ylabel('Végétation');
title('diagramme de V en fonction de P');
grid on
hold off
