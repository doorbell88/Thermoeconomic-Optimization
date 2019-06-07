function [T,p,G,Z,efficiency,cost] = newParameter(x);

%-------------------------------------------------------------------------------
% DECISION VARIABLES
beta_c = x(1);  % Compressor compression ratio
eta_c  = x(2);  % Compressor isentropic efficiency
T(3)   = x(3);  % Turbine inlet temperature
eta_t  = x(4);  % Turbine isentropic efficiency
%-------------------------------------------------------------------------------

% Givens
W_EL_NET = 20000; %[kW]
eff_gen  = 0.98;  % Generator efficiency

%-------------------------- Thermodynamic Parameters ---------------------------
% Ambient conditions
T0 = 25 + 273.15; %[K]
p0 = 1; %[bar]

cp_a = 1.004;   % heat capacity of air
cp_g = 1.17;    % heat capacity of exhaust gas

gamma_a = 1.4;  % adiabatic coefficient of air
gamma_g = 1.33; % adiabatic coefficient of exhaust gas

LHV = 50000;    % LHV of CH4 [kJ/kg]
rb  = 0.95;     % Combustion chamber Pout/Pin


%===============================================================================
%--------------------------------- Compressor ----------------------------------
% Temperature at outlet of compressor
k_a = (gamma_a - 1) / gamma_a;
T(2) = T0*(1 + (((beta_c^k_a)-1) / eta_c));

% Pressure at outlet of compressor
p(2) = p0 * beta_c;


%----------------------------- Combustion Chamber ------------------------------
T(1) = T0;
p(1) = p0; %p(2);
p(3) = p(2)*rb;
T(3) = x(3);    % (DECISION VARIABLE)


%---------------------------------- Turbine ------------------------------------
% pressure at turbine outlet
p(4) = p0;          % pressure at turbine outlet
rt = p(3)/p(4);     % expansion ratio

% Temperature at outlet of turbine
k_g = ((gamma_g-1)/gamma_g);
T(4) = T(3) * (1 - (eta_t*(1- ((1/rt)^k_g) )) );


%===============================================================================
%--------------------------------- Mass Flows ----------------------------------

% Fuel to air ratio
f = (cp_g*(T(3)-T0) - cp_a*(T(2)-T0)) / (LHV - cp_g*(T(3)-T0));

% Mass flowrate of air
Ga = (W_EL_NET/eff_gen) / ( cp_g*(1+f)*(T(3)-T(4)) - cp_a*(T(2)-T0) );

% Mass flowrate of fuel
Gf = f*Ga;

% Mass flowrate of exhaust gas
Gg = Ga + Gf;

% Mass flowrate of CO
% https://www3.epa.gov/ttnchie1/ap42/ch03/final/c03s01.pdf
kgCO_kgCH4 = 1.77e-3;   %[kgCO/kgCH4]
Gco = kgCO_kgCH4 * Gf;  %[kgCO]

% Compile mass flow vector
G(1) = Gf;
G(2) = Ga;
G(3) = Gg;
G(4) = G(3);
G(5) = Gco;


%--------------------------------- ECONOMICS -----------------------------------
% NOMINAL OPERATION
phi = 1.05;
tau = 3600; %[s/hr]
N   = 8000; %[hrs/yr]
CRF = 0.25;
cf  = 5.00e-6; %[$/kJ]

%---------------------------------------
% EMISSIONS
z_co   = 10;    %[$/kgCO]
fp_co  = 6.25;  %[-]
M_co   = 28.01; %[kgCO/kmolCO]
M_ch4  = 16.04; %[kgCH4/kmolCH4]

%---------------------------------------
% ECONOMIC PARAMETERS
c11 = 39.5; %[$/(kg/s)]
c12 = 0.9;

c21 = 25.6; %[$/(kg/s)]
c22 = 0.995;
c23 = 0.018; %[1/K]
c24 = 26.4;

c31 = 266.3; %[$/(kg/s)]
c32 = 0.92;
c33 = 0.036; %[1/K]
c34 = 54.4;


%------------------------------- Cost Functions --------------------------------
% Compressor investment cost
Z(1) = ((c11*Ga)/(c12-eta_c)) * beta_c*log(beta_c);

% Combustor investment cost
Z(2) = ((c21*Ga)/(c22-rb)) * (1 + exp(c23*T(3) - c24));

% Gas Turbine investment cost
Z(3) = ((c31*Gg)/(c32-eta_t)) * log(rt) * (1 + exp(c33*T(3) - c34));


%---------------------------------------
% Total investment cost [$]
Z(5) = Z(1)+Z(2)+Z(3);

%---------------------------------------
% Operational cost [$/s]
Z(6) = (Z(5)*CRF*phi) / (tau*N);

% Fuel cost [$/s]
Z(7) = cf*Gf*LHV;

% Emission Cost [$/s]
Z(8) = fp_co * z_co * Gco;


%--------------------------------- Efficiency ----------------------------------
efficiency = W_EL_NET / (Gf * LHV);

%---------------------------------------
% FINAL OBJECTIVE FUNCTION
cost_per_second = Z(6)+Z(7)+Z(8);           %[$/s]
cost_per_year   = cost_per_second*tau*N;    %[$/yr]

cost = cost_per_year;
end
