clear all
% Initialization
[~,sys_eq,x0,u0] = s6stgPSS;
% Approximation
[~,A_l, B_l, A_lq, B_lq, N_lq, A_lqc, B_lqc, N_lqc, ids_lq, ids_lqc] = smibaprox(sys_eq,x0,u0);
% Reducing
%[A_lq, B_lq, N_lq, ids_lq] = smibreduce(A_lq,B_lq,N_lq,length(x0));
%[A_lqc, B_lqc, N_lqc, ids_lqc] = smibreduce(A_lqc,B_lqc,N_lqc,length(x0));
% Dimensions of the systems
sl = size(A_l,1);
slq = size(A_lq,1);
%% Simulation
% (a) settings:
t_init = 0; % initial time of simulation
t_fin = 15; % final time of simulation
dt = 1e-4; % time step
% Perturbation settings:
nxp = [2]; % number of perturbed state
vxp = 0.0; % state perturbation magnitude
nup = 1; % number of perturbed input
vup = 1.6; % input perturbation magnitude
tu = [  01.00   01.10;  % time of first input perturbation
        00.00   00.00   % time of second input perturbation
                    ];
t = [t_init t_fin];
x_per = zeros(size(x0));
u_per = zeros(size(u0));
x_per(nxp) = vxp;
u_per(nup) = vup;
% (b) time domain response:
[t,y_l,y_lq,y_lqc,y_nl] = smibsim (t,dt,A_l,B_l,A_lq,B_lq,N_lq,A_lqc,...
    B_lqc,N_lqc,sys_eq,x0,u0,x_per,u_per,tu,ids_lq,ids_lqc);
%% (c) plots (compare approximations)
for nst = 1:sl % number of compared state
figure('WindowStyle', 'docked', 'ToolBar', 'none')
% plotting:
plot(t,y_nl(:,nst)-x0(nst).','b','LineWidth',1.0)
hold on
plot(t,y_lqc(:,nst),'c','LineWidth',1.0)
hold on
plot(t,y_lq(:,nst),'r','LineWidth',1.0)
hold on
plot(t,y_l(:,nst),'g','LineWidth',1.0)
hold off
% scale
xlim([1 5])
% decor
legend('nl','lqc','lq','lin')
%plotbrowser('on')
%propertyeditor('on')
grid on
grid minor
xlabel ('Время, с')
title (['\DeltaT_m = ',num2str(vup)])
%title(['ns=',num2str(ns),'; nm=',num2str(nm),'; x=[',num2str(x_per(1)), '  ',...
%    num2str(x_per(2)), '  ',num2str(x_per(3)),']'])
end
%% Errors of approximations
% Simulation settings:
t_init = 0; % initial time of simulation
t_fin = 5; % final time of simulation
dt = 1e-4; % time step
% Perturbation settings:
nxp = [2]; % number of perturbed state
vxp = 0.0; % state perturbation magnitude
nup = 2; % number of perturbed input
vup = [0.01:0.01:1.7]; % input perturbation magnitude
tu = [  01.00   01.10;  % time of first input pertirbation
        00.00   00.00   % time of second input pertirbation
                    ];
x_per = zeros(size(x0));
u_per = zeros(size(u0,1),length(vup));
x_per(nxp) = vxp;
u_per(nup,:) = vup;
nbytes = fprintf('Time Domain Simulation: 00%% \n');
tot = length(vup);
cur = 0;
x0t = x0.';
L2_nl_lin = [];
L2_nl_lq = [];
L2_nl_lqc = [];
% Time Domain Simulation:
tic
for k = 1:length(vup)
    cur = cur + 1;
    t = [t_init t_fin];
    [t,y_l,y_lq,y_lqc,y_nl]=smibsim (t,dt,A_l,B_l,A_lq,B_lq,N_lq,A_lqc,...
        B_lqc,N_lqc,sys_eq,x0,u0,x_per,u_per(:,k),tu,ids_lq,ids_lqc);
    L2_nl_lin(k,:) = trapz(t,(y_nl-x0t-y_l).^2);
    L2_nl_lq(k,:) = trapz(t,(y_nl-x0t-y_lq(:,1:sl)).^2);
    L2_nl_lqc(k,:) = trapz(t,(y_nl-x0t-y_lqc(:,1:sl)).^2);
    %L2_lq_lin(k,:) = trapz(t,(y_lq(:,1:sl)-y_l).^2);
    step = cur / tot;
    fprintf(repmat('\b',1,nbytes))
    nbytes = fprintf('Time Domain Simulation: %d%% \n',round(step*100));
end
toc
% Plot of errors
for nst = 2
    figure()
    plot(vup,sqrt(L2_nl_lin(:,nst)),vup,sqrt(L2_nl_lq(:,nst)),vup,...
        sqrt(L2_nl_lqc(:,nst)),'LineWidth',1.5)
    grid on
    grid minor
    legend('лин.', 'квадр.', 'куб.')
    %xlabel ('\DeltaE_{fd}')
    ylabel('Ошибка по \delta')
end