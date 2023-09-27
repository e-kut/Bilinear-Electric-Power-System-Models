function [sys_init,sys_eq,x0,u0] = s2stg(varargin)
%----------------------------- Initialization -----------------------------
%                       Example 12.2 p.732, Kundur

% Per Unit System: 2220 MVA, 24 kV
P = 0.9;            % active power, [pu]
Q = 0.3;            % reactive power, [pu]
E_tmod = 1.0;       % terminal voltage magnitude, [pu]
E_tdph = 36;        % terminal voltage phase in degrees
E_Bmod = 0.995;     % infinite bus voltage magnitude, [pu]
E_Bdph = 0;         % infinite bus voltage phase in degrees
X_d_ht = 0.3;       % transient reactance of the d-axis stator winding
X_T = 0.95;         % total system reactance, [pu]
H = 3.5;            % inertia moment, [MW*s/MVA]
sys_freq = 60;      % system frequency, [Hz]
K_D = 10;           % damping coefficient
%------------------------- End of Initialization --------------------------


%------------------- Calculating the Operating Point ----------------------
% stedy-state value of the rotor speed:
domega_0 = 0;
% complex terminal voltage with E_t as a reference phasor (angle E_t = 0):
E_tc = E_tmod * exp(1i*deg2rad(0));
% generator stator current:
I_tc = conj(P + 1i*Q) / conj(E_tc);
% the voltage behind the transient reactance:
E_htc = E_tc + 1i*X_d_ht * I_tc;
% the angle by wich E_htc leads E_Bc (stedy-state value of rotor angle):
delta_0d = rad2deg(angle(E_htc)) + E_tdph; % in degrees
delta_0r = deg2rad(delta_0d); % in radians 
% state vector:
x0 = [domega_0; delta_0r];
% inpiut variable:
T_m0 = abs(E_htc) * E_Bmod * sin(delta_0r) / X_T;
% inpiut vector:
u0 = T_m0;
%--------------- End of Calculation of Operating Point --------------------


%------------------- Power System Dynamics Equations ----------------------
% Constructing a function for electrical torque
Te = @(x) abs(E_htc) * E_Bmod * sin(x(2)) / X_T;
% Construction of generator dynamic equations
dyn_eq{1,1} =  @(x) -K_D / (2 * H) * x(1) - Te(x) / (2 * H);
dyn_eq{2,1} = @(x) 2*pi*sys_freq * x(1);
% Construction of input equations
inp_eq = @(u) u / (2*H);
% Construction of the overall EPS equations
sys_eq {1,1} = @(x,u) dyn_eq{1,1}(x) + inp_eq(u);
sys_eq {2,1} = @(x,u) dyn_eq{2,1}(x);
%----------------- End of Power System Dynamics Equations -----------------


%----------------------- Construct Output Structure -----------------------
sys_init.sys_eq = sys_eq;
sys_init.x0 = x0;
sys_init.u0 = u0;
%------------------- End of Construct Output Structure --------------------
end
