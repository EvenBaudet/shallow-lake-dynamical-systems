c = constants(); 

Veq = 0.301;
Peq = 0.103;

% Jacobian, 2D-system

function dfvdV = dfvdV(V, Pw, c)
    term = (c.gammae*Pw*c.hv/c.he)^c.pe;
    dfvdV = c.rv - 2*V*term*(1/(c.hv+V))^c.pe + (V^2)*term*c.pe*(1/(c.hv+ ...
        V)^(c.pe+1));
end

function dfvdP = dfvdP(V, Pw, c)
    term = (c.gammae*Pw*c.hv/c.he)^c.pe;
    dfvdP = -(V^2) *(term*(1/(c.hv+V))^c.pe)*c.pe*Pw^(c.pe-1);
  
end

% Calcul 

dfvdV_eq = dfvdV(Veq, Peq, c);
dfvdP_eq = dfvdP(Veq,Peq, c);
dfpdP_eq = 1/c.s;
dfpdV_eq = 0;

detJ = dfvdV_eq*dfpdP_eq;

fprintf('J11 : %f\n', dfvdV_eq);
fprintf('J22 : %f\n', dfpdP_eq);
fprintf('J12 : %f\n',dfvdP_eq);
fprintf('detJ : %f\n',detJ);
