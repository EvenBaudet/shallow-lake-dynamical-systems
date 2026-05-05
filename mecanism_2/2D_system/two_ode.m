function dXdt = two_ode(t,X)

c = constants();

V  = X(1);
Pw = X(2);

% --- SOD equilibrium

sod_eq = c.kv*V/c.l;

% --- Turbidity

e = c.gammae*Pw*(c.hv/(c.hv+V));

% --- Hill

Hill = (sod_eq^c.psod) / (sod_eq^c.psod + c.hsod^c.psod);

% --- Differential equations

dVdt = c.rv*V*(1-V*(((c.he^c.pe)+(e^c.pe))/c.he^c.pe));
dPwdt = (c.pwin-Pw)/c.s + c.g*Hill;

dXdt = [dVdt; dPwdt];

end
