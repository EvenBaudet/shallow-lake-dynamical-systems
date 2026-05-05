function nullclines_and_eq()
    
    c = constants();
    
    % Data
    V_vec = linspace(0.001, 0.999, 500);
    
    % P nullcline Pw = Pwin + τ*g*Hill(V)
    Pw_V_null = (c.he / c.gammae)*((1 - V_vec)./V_vec).^(1/c.pe).*((c.hv + ...
        V_vec)/c.hv);
    
    % SOD term
    SOD_vec = (c.kv / c.l) * V_vec;

    Hill_vec = (SOD_vec.^c.psod) ./ (SOD_vec.^c.psod + c.hsod^c.psod);

    % Calcul Pw
    Pw_P_null = c.pwin + c.s * c.g * Hill_vec;
    
    
    % V tq Pw_V_null(V) - Pw_P_null(V) = 0
    diff_func = @(V) compute_Pw_V(V, c) - compute_Pw_P(V, c);
    initial_guesses = [0.1, 0.3, 0.5, 0.7, 0.9];
    equilibria_V = [];
    equilibria_Pw = [];
    
    for guess = initial_guesses
        try
            V_eq = fzero(diff_func, guess);
            
            Pw_eq = compute_Pw_V(V_eq, c);
           
            if V_eq > 0 && V_eq < 1 && Pw_eq > 0
                is_new = true;
                for i = 1:length(equilibria_V)
                    if abs(V_eq - equilibria_V(i)) < 1e-4 && abs(Pw_eq - equilibria_Pw(i)) < 1e-4
                        is_new = false;
                        break;
                    end
                end
                
                if is_new
                    equilibria_V = [equilibria_V, V_eq];
                    equilibria_Pw = [equilibria_Pw, Pw_eq];
                end
            end
        catch
            continue;
        end
    end
    
 
    % Plot
    
    figure;
    hold on;
    grid on;
    box on;
    
    % V nullcline
    plot(Pw_V_null, V_vec, 'g.', 'LineWidth', 3, 'DisplayName', ...
        'Nullcline V (dV/dt=0)');
    
    % Pw nullcline
    plot(Pw_P_null, V_vec, 'r.', 'LineWidth', 3, 'DisplayName', ...
        'Nullcline Pw (dPw/dt=0)');
    
    % Eq
    if ~isempty(equilibria_V)
        scatter(equilibria_Pw, equilibria_V, 150, 'k', 'filled', ...
                'DisplayName', 'Points d''équilibre', 'LineWidth', 2);
        
        for i = 1:length(equilibria_V)
            text(equilibria_Pw(i) + 0.002, equilibria_V(i) + 0.02, ...
                 sprintf('E%d', i), 'FontSize', 12, 'FontWeight', 'bold');
        end
    end
    
    % Pw(V)
    try
        [t, X, V_sim, Pw_sim] = run_simulation();
        
        % We only keep the data after 5 years
        t_years = t/365;
        idx = t_years > 5;
        plot(Pw_sim(idx), V_sim(idx), 'b-', 'LineWidth', 1.5, ...
             'DisplayName', 'Trajectoire');
    catch
        warning('La simulation a échoué, tracé sans trajectoire');
    end
    
    
    xlabel('P_w (g m^{-3})', 'FontSize', 14);
    ylabel('V', 'FontSize', 14);
    title('PNullclines of V and P with the equilibrium point', 'FontSize', ...
        16);
    legend('Location', 'best', 'FontSize', 12);
    
    xlim([0, max([Pw_V_null, Pw_P_null]) * 0.9]);
    ylim([0, 1]);
 
    grid on;
    grid minor;
    
    if ~isempty(equilibria_V)
        subtitle_str = sprintf('%d points d''équilibre trouvés', length(equilibria_V));
        for i = 1:length(equilibria_V)
            subtitle_str = sprintf('%s\nE%d: V=%.3f, Pw=%.4f', subtitle_str, i, equilibria_V(i), equilibria_Pw(i));
        end
        subtitle(subtitle_str, 'FontSize', 10);
    end
    
    hold off;
    
    % equilibres
    
    disp('=== POINTS D''ÉQUILIBRE ===');
    for i = 1:length(equilibria_V)
        disp(sprintf('E%d: V = %.4f, Pw = %.5f g/m³', i, equilibria_V(i), equilibria_Pw(i)));
    end
end

% ========== TOOLS ==========

function Pw = compute_Pw_V(V, c)
    % Calcul Pw on the nullcline
    if V == 0
        Pw = Inf;  
    else
        Pw = (c.he / c.gammae) * ((1 - V)/V)^(1/c.pe) * ((c.hv + V)/c.hv);
    end
end

function Pw = compute_Pw_P(V, c)
    % Calcule Pw on the Pw nullcine for a V
    SOD = (c.kv / c.l) * V;
    Hill = (SOD^c.psod) / (SOD^c.psod + c.hsod^c.psod);
    Pw = c.pwin + c.s * c.g * Hill;
end

function [t, X, V, Pw] = run_simulation()
    t0 = 0;
    tf = 370*25;  
    tspan = [t0 tf];
    
    % Initial values
    V0 = 0.2;
    Pw0 = 0.05;
    SOD0 = 1;
    X0 = [V0; Pw0; SOD0];
    
    % Resolution
    [t, X] = ode15s(@three_ode, tspan, X0);
    
    % Extraction variables
    V = X(:,1);
    Pw = X(:,2);
end