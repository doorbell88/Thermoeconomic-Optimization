clear all; clc;

%-------------------------------------------------------------------------------
% DECISION VARIABLES
% beta_c = x(1);  % Compressor compression ratio
% eta_c  = x(2);  % Compressor isentropic efficiency
% T(3)   = x(3);  % Combustion chamber exit temperature
% eta_t  = x(4);  % Turbine isentropic efficiency
%-------------------------------------------------------------------------------
%    (1)        (2)         (3)         (4)  
%    beta_c     eta_c       T3          eta_t
%-------------------------------------------------
x0 = [17,       0.85,       1400,       0.85    ];
lb = [10,       0.8,        1000,       0.8     ];
ub = [25,       0.89,       1800,       0.91    ];


%-------------------------------------------------------------------------------

% Get options for fmincon from optimset
options = optimset('algorithm', 'interior-point');

% Optimization algorithm
x = fmincon('objF', x0, [],[],[],[], lb, ub, 'constraints', options);

% give decision variables
[T,p,G,Z,efficiency,cost] = newParameter(x);

% print out final result
format compact
disp("====================================================")
disp("VECTORS")
disp("-------")
disp(x)
T
p
G
Z

disp("====================================================")
disp("DECISION VARIABLES")
disp("------------------")
disp("beta_c, Compressor compression ratio")
disp(x(1))

disp("eta_c, Compressor isentropic efficiency")
disp(x(2))

disp("T(3), Combustion chamber exit temperature")
disp(x(3))

disp("eta_t, Turbine isentropic efficiency")
disp(x(4))

% disp("====================================================")
% disp("TEMPERATURES")
% disp("------------")
% disp("T(1)")
% disp(T(1))
% 
% disp("T(2)")
% disp(T(2))
% 
% disp("T(3)")
% disp(T(3))
% 
% disp("T(4)")
% disp(T(4))
% 
disp("====================================================")
disp("MASS FLOWRATES [kg/s]")
disp("---------------------")
disp("Air [kg/s]")
disp(G(2))
disp("Fuel [kg/s]")
disp(G(1))
disp("Gas [kg/s]")
disp(G(3))
disp("Carbon Monoxide (CO) [kg/s]")
disp(G(5))


disp("====================================================")
disp("THERMODYNAMIC CYCLE EFFICIENCY")
disp(efficiency)


W_EL_NET = 20000; %[kW]
tau = 3600; %[s/hr]
N   = 8000; %[hrs/yr]
disp("====================================================")
disp("INVESTMENT AND OPERATIONAL COSTS")
disp("--------------------------------")
disp("Compressor investment cost [$]")
disp(Z(1))

disp("Combustor investment cost [$]")
disp(Z(2))

disp("Gas Turbine investment cost [$]")
disp(Z(3))

disp("______________________________")
disp("Total investment cost [$]")
disp(Z(5))

disp("..............................")
disp("Recovery cost [$/s], [$/yr]")
disp(Z(6))
disp(Z(6)*tau*N)

disp("Fuel cost [$/s], [$/yr]")
disp(Z(7))
disp(Z(7)*tau*N)

disp("Emission Cost [$/s], [$/yr]")
disp(Z(8))
disp(Z(8)*tau*N)


disp("====================================================")
disp("FINAL LEVELIZED COST [USD/yr]")
disp(cost)
disp("Levelized cost of electricity [$/kWh]")
LCOE = cost/(W_EL_NET*N);
disp(LCOE)
