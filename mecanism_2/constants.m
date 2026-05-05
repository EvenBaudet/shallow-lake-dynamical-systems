function c = constants()

% --- Vegetation

c.rv = 0.07;    % Maximum growth rate of vegetation 
c.kv = 0.05;    % Maximum attribution of vegetation to SOD
c.hv = 0.2;     %Vegetation cover with half of the maximal 
                %effect on water clarity

% --- Phosphore

c.pwin = 0.05;  % Phosphorus concentration incoming water
c.s = 475;      % Retention time (day) : s augmente, période augmente

% --- SOD

c.l = 0.01;     % Mineralization rate of sediment oxygen demand (day) 
c.g = 0.002;    % Maximal effect of SOD on phosphorus
c.hsod = 2;     % Critical sediment oxygen demand

% --- Turbidity

c.he = 1;
c.gammae = 30;
c.pe = 5;
c.psod = 10;

end
