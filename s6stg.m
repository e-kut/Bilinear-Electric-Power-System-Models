function [sys_init,sys_eq,x0,u0] = s6stg(varargin)
%----------------------------- Initialization -----------------------------
%               Example 12.3 p.752, Example 12.5, p.789 Kundur

% Per Unit System: 2220 MVA, 24 kV
sys_freq = 60;  % system frequency, [Hz]
X_d = 1.81;     % synchronous d-axis stator winding reactance, [p.u]
X_d_ht = 0.3;   % transient d-axis stator winding reactance, [p.u]
X_d_2ht = 0.23; % subtransient d-axis stator winding reactance, [p.u]
X_q = 1.76;     % synchronous q-axis stator winding reactance, [p.u]
X_q_ht = 0.65;  % transient q-axis stator winding reactance, [p.u]
X_q_2ht = 0.25; % subtransient q-axis stator winding reactance, [p.u]
T_d0_ht = 8;    % transient d-axis open-circuit time constant, [sec]
T_d0_2ht = 0.03;% subtransient d-axis open-circuit time constant, [sec]
T_q0_ht = 1;    % transient q-axis open-circuit time constant, [sec]
T_q0_2ht = 0.07;% subtransient q-axis open-circuit time constant, [sec]
R_a = 0.003;    % stator winding resistance, [p.u]
X_l = 0.16;     % stator winding leakage inductance, [p.u]
H = 3.5;        % inertia moment, [MW*s/MVA]
K_D = 0;        % damping coefficient
R_E = 0;        % line resistance, [p.u]
X_E = 0.65;     % line reactance, [p.u]
A_Sat = 0.031;  % saturation characteristic
B_Sat = 6.93;   % saturation characteristic
psi_TI = 0.8;   % saturation characteristic
P = 0.9;        % active power, [p.u]
Q = 0.3;        % reactive power, [p.u]
E_t = 1.0;      % terminal voltage magnitude, [p.u]
%------------------------- End of Initialization --------------------------

%------------------- Calculating the Operating Point ----------------------
omega_0 = 2*pi*sys_freq;
% Saturation mode identifying
if ~isempty(varargin)
    sat_mod = varargin{1};
else
    sat_mod = string;
end
% Stator winding inductances (in per unit system X_index = L_index):
L_d = X_d;
L_d_ht = X_d_ht;
L_d_2ht = X_d_2ht;
L_q = X_q;
L_q_ht = X_q_ht;
L_q_2ht = X_q_2ht;
L_l = X_l;
% Unsaturated mutual inductances(Example 4.1 (p.154)):
L_ad = L_d - L_l;
L_aq = L_q - L_l;
% Field winding inductance & resistance(Example 4.1,p.154,Table 4.1,p.150):
L_fd = (L_d - L_l) * (L_d_ht - L_l) / (L_d - L_d_ht);
R_fd = (L_d - L_l)^2 / ((L_d - L_d_ht) * T_d0_ht * omega_0);
% Damping windings inductances (Example 4.1 p.154, Table 4.1, p.150):
L_1q = (L_q - L_l) * (L_q_ht - L_l) / (L_q - L_q_ht);
L_1d = -(L_d_2ht - L_l) * (L_d_ht - L_l) / (L_d_2ht - L_d_ht);
L_2q = -(L_l - L_q_2ht) * (L_l - L_q_ht) / (L_q_2ht - L_q_ht);
% Damping windings resistances (Example 4.1 p.154, Table 4.1, p.150):
R_1d = -(L_l - L_d_ht)^2 / ((L_d_2ht - L_d_ht) * T_d0_2ht * omega_0);
R_1q = (L_l - L_q)^2 / ((L_q - L_q_ht) * T_q0_ht * omega_0);
R_2q = -(L_l - L_q_ht)^2 / ((L_q_2ht - L_q_ht) * T_q0_2ht * omega_0);
% Total saturation factors K_sd, K_sq (Section 3.8):
I_t_tild = (P + Q * 1i)' / E_t'; % terminal current, p.733
if sat_mod == "nosat"
    K_sd = 1;
    K_sq = 1;
else
    % air-gap voltage:
    E_a_tild = E_t + (R_a + X_l * 1i) * I_t_tild; % eq.3.194, p.116
    % air-gap flux linkage:
    psi_at = abs(E_a_tild);
    % deviation psi_at from psi_at0:
    if psi_at <= psi_TI
        psi_I = 0; % eq.3.188, p.116
    else
        psi_I = A_Sat * exp(B_Sat * (psi_at-psi_TI)); % eq.3.189, p.116
    end
    % saturation factors:
    K_sd = psi_at / (psi_at + psi_I); % eq.3.187, p.116
    K_sq = K_sd;
end
% Saturated values of the mutual inductances without field and damping 
% windings:
L_ads = K_sd * L_ad;
L_aqs = K_sq * L_aq;
% Saturated values of synchronous reactances & inductances (p.746):
X_ds = L_ads + L_l;
X_qs = L_aqs + L_l;
L_ds = X_ds;
L_qs = X_qs;
% Steady-state value of the terminal current & power factor angle:
I_t = sqrt(P^2 + Q^2) / E_t; % (p.101)
f_r = acos(P / (E_t * I_t)); % in radians (p.101)
f_d = f_r * 180/pi; % in degrees
% Steady-state value of the internal angle (p.746):
num = I_t * X_qs * cos((f_r)) - I_t * R_a * sin((f_r));
denom = E_t + I_t * R_a * cos((f_r)) + I_t * X_qs * sin((f_r));
delta_ir = atan(num / denom); % in radians in the interval [-pi/2, pi/2]
% delta_ir = delta_ir + 0.5 * (1 - sign(denom)) * pi + 0.5 * (1 + ...
%     sign(denom)) * (1 - sign(num)) * pi; % in radians in the interval [0, 2*pi]
delta_id = delta_ir * 180 / pi; % in degrees
% Steady-state value of the rotor (stator?) voltages & currents (p.746):
e_d0 = E_t * sin (delta_ir);
e_q0 = E_t * cos (delta_ir);
i_d0 = I_t * sin (delta_ir + f_r); % only if f is lagging power factor angle
i_q0 = I_t * cos (delta_ir + f_r); % only if f is lagging power factor angle
% Steady-state value of the bus voltage (p.746):
E_Bd0 = e_d0 - R_E * i_d0 + X_E * i_q0;
E_Bq0 = e_q0 - R_E * i_q0 - X_E * i_d0;
E_B = sqrt(E_Bd0^2 + E_Bq0^2);
% Steady-state value of the rotor angle (p.746):
delta_0r = atan(E_Bd0 / E_Bq0); % in radians in the interval [-pi/2, pi/2]
delta_0r = delta_0r + 0.5 * (1 - sign(E_Bq0)) * pi + 0.5 * (1 + sign(E_Bq0))...
    * (1 - sign(E_Bd0)) * pi; % in radians in the interval [0, 2*pi]
delta_0 = delta_0r * 180 / pi; % in degrees
% Steady-state value of the field winding current & voltage (p.746):
i_fd0 = (e_q0 + R_a * i_q0 + L_ds * i_d0) / L_ads;
E_fd0 = L_ad * i_fd0;
% Steady-state value of the mutual flux linkage (p.746):
psi_ad0 = L_ads * (-i_d0 + i_fd0);
psi_aq0 = - L_aqs * i_q0;
% Steady-state value of the field & damping windings flux linkages (p.104):
psi_fd0 = (L_ads + L_fd) * i_fd0 - L_ads * i_d0;
psi_1d0 = L_ads * (i_fd0 - i_d0);
psi_1q0 = -L_aqs * i_q0;
psi_2q0 = -L_aqs * i_q0;
% Steady-state value of the mechanical torque:
T_m0 = -(psi_aq0 * i_d0 - psi_ad0 * i_q0);
% Operating point:
x0 = [0; delta_0r; psi_fd0; psi_1d0; psi_1q0; psi_2q0]; % states
u0 = [T_m0; E_fd0]; % inputs
% Incremental saturation (p.745)
if sat_mod == "incsat"
    psi_at0 = psi_at;
    % incremental saturation factor:
    K_sd_incr = 1 / (1 + B_Sat * A_Sat * exp(B_Sat * (psi_at0 - psi_TI)));
    K_sq_incr = K_sd_incr;
    % incremental saturated values of the mutual inductances:
    L_ads = K_sd_incr * L_ad;
    L_aqs = K_sq_incr * L_aq;
end
% Saturated values of the mutual inductances with field and damping
% windings:
L_ads_2ht = 1 / (1 / L_ads + 1 / L_fd + 1 / L_1d);
L_aqs_2ht = 1 / (1 / L_aqs + 1 / L_1q + 1 / L_2q);
%--------------- End of Calculation of Operating Point --------------------

%------------------- Power System Dynamics Equations ----------------------
% Auxiliary coefficients:
X_Td = X_E + L_ads_2ht + L_l;
X_Tq = X_E + L_aqs_2ht + L_l;
R_T = R_a + R_E;
K_XR = R_T^2 + X_Td * X_Tq;
% Auxiliary functions:
% (a) d- & q-axis stator currents:
i_d = @(x) X_Tq / K_XR * (L_ads_2ht * (x(3) / L_fd + x(4) / L_1d) - ...
    E_B * cos(x(2))) - R_T / K_XR * (L_aqs_2ht * (x(5) / L_1q + ...
    x(6) / L_2q) + E_B * sin(x(2)));
i_q = @(x) R_T / K_XR * (L_ads_2ht * (x(3) / L_fd + x(4) / L_1d) - ...
    E_B * cos(x(2))) + X_Td / K_XR * (L_aqs_2ht * (x(5) / L_1q + ...
    x(6) / L_2q) + E_B * sin(x(2)));
% (b) d & q-axis mutual flux linkage
psi_ad = @(x) L_ads_2ht * (x(3) / L_fd + x(4) / L_1d - i_d(x));
psi_aq = @(x) L_aqs_2ht * (x(5) / L_1q + x(6) / L_2q - i_q(x));
% (c) field & damping winding currents:
i_fd = @(x) (x(3) - psi_ad(x)) / L_fd;
i_1d = @(x) (x(4) - psi_ad(x)) / L_1d;
i_1q = @(x) (x(5) - psi_aq(x)) / L_1q;
i_2q = @(x) (x(6) - psi_aq(x)) / L_2q;
% States dynamic equations:
dyn_eq{1,1} = @(x) 1 / (2*H) * (-K_D * x(1) - (psi_ad(x) * i_q(x) -...
    psi_aq(x) * i_d(x)));
dyn_eq{2,1} = @(x) omega_0 * x(1);
dyn_eq{3,1} = @(x) -omega_0 * R_fd * i_fd(x);
dyn_eq{4,1} = @(x) -omega_0 * R_1d * i_1d(x);
dyn_eq{5,1} = @(x) -omega_0 * R_1q * i_1q(x);
dyn_eq{6,1} = @(x) -omega_0 * R_2q * i_2q(x);
% Inputs (controls):
inp_eq {1,1} = @(u) u(1) / (2 * H);
inp_eq {2,1} = @(u) u(2) * omega_0 * R_fd / L_ad;
% Overall system equations:
sys_eq{1,1} = @(x,u) dyn_eq{1,1}(x) + inp_eq{1,1}(u);
sys_eq{2,1} = @(x,u) dyn_eq{2,1}(x);
sys_eq{3,1} = @(x,u) dyn_eq{3,1}(x) + inp_eq{2,1}(u);
sys_eq{4,1} = @(x,u) dyn_eq{4,1}(x);
sys_eq{5,1} = @(x,u) dyn_eq{5,1}(x);
sys_eq{6,1} = @(x,u) dyn_eq{6,1}(x);
%----------------- End of Power System Dynamics Equations -----------------
sys_init.sys_eq = sys_eq;
sys_init.x0 = x0;
sys_init.u0 = u0;
%------------------------------ Phasor Plot -------------------------------
% plot([0,e_d0],[0,e_q0])
% hold on
% plot([0,i_d0],[0,i_q0])
% hold on
% plot([0,E_Bd0],[0,E_Bq0])
% hold on
% plot([0,P],[0,Q])
% hold off
% xlm = 2;
% ylm = xlm;
% xlim([-xlm xlm])
% ylim([-ylm ylm])
% xlabel ('d-axis')
% ylabel ('q-axis')
% grid on
% grid minor
% legend ('E_t', 'I_t', 'E_B', 'PQ')
%--------------------------- End of Phasor Plot ---------------------------
end