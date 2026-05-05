function dXdt = three_ode(t,X)

c = constants();

V  = X(1);
Pw = X(2);
SOD = X(3);

% --- Turbidity

e = c.gammae*Pw*(c.hv/(c.hv+V));

% --- Hill

Hill = (SOD^c.psod) / (SOD^c.psod + c.hsod^c.psod);

% --- Differential equations

dVdt = c.rv*V*(1-V*(((c.he^c.pe)+(e^c.pe))/c.he^c.pe));
dPwdt = (c.pwin-Pw)/c.s + c.g*Hill;
dSODdt = c.kv*V - c.l*SOD;

dXdt = [dVdt; dPwdt; dSODdt];

end