function [t,y_l,y_lq,y_lqc,y_nl] = smibsim (t,dt,A_l,B_l,A_lq,B_lq,N_lq,A_lqc,...
    B_lqc,N_lqc,sys_eq,x0,u0,x_per,u_per,tu,ids_lq,ids_lqc,varargin)
%------------------------------ Description -------------------------------
%
% This function performs time domain simulation of the linear, 
% linear-quadratic and nonlinear SMIB EPS
%
% INPUTS:   t       -   (row-vector) time limits of the simulation in the
%                       following format t = [start-time  end-time]
%           dt      -   (scalar) time step
%           u0      -   (column-vector) initial steady-state values of the
%                       input variables
%           A_l     -   dynamics matrix of the linear system
%           B_l     -   input matrix of the linear system
%           A_lq    -   dynamics matrix of the linear-quadratic system
%           B_lq    -   input matrix of the linear-quadratic system
%           N       -   (structure) a set of matrices of bilinear terms
%           sys_eq  -   (cell array) set of the system equations in handle 
%                       function form
%           x0      -   (column-vector) initial steady-state values of the
%                       state variables
%           u0      -   (scalar or column-vector) initial steady-state 
%                       values of the input variables
%           x_per   -   (column-vector) deviation from x0
%           u_per   -   (scalar or column-vector) deviation from u0
%           tu      -   start and end time of submission perturbation 
%                       through the system input (u_per). For a scalar 
%                       u_per tu is the row-vector with the following 
%                       format tu = [start-time  end-time]. For a vector 
%                       u_per of length n tu is the matrix with following 
%                       format:
%                       tu = [start-time-for-u_per(1) end-time-for-u_per(1);
%                            start-time-for-u_per(2) end-time-for-u_per(2);
%                            .............................................;
%                            start-time-for-u_per(n) end-time-for-u_per(n)]
%           ids     -   numbers of saved states in the original vector of 
%                       state variables after reducion. If the reduction 
%                       was not performed then ids must be equal to vector
%                       [1 : n+n^2] where n is the number of states in the 
%                       linear system
%           varargin =  1 - only linear system simulation
%                       2 - only linear-quadratic system simulation
%                       3 - only linear-quadratic-cubic system simulation
%                       4 - only nonlinear system simulation
%                       0 - simulation of the all systems
%
% OUTPUTS:  t       -   time
%           y_l     -   dynamic of linear system states
%           y_lq    -   dynamic of linear-quadratic system states
%           y_lqc   -   dynamic of linear-quadratic-cubic system states
%           y_nl    -   dynamic of nonlinear system states
%
% Linear System: dx/dt = A_l * x + B_l * u
% Linear-Quadratic System: dx/dt = A_lq * x + sum_i(Ni*x*ui) + B_lq * u
%--------------------------- End of Description ---------------------------

%----------------------------- Initialization -----------------------------
t = t(1):dt:t(2); % time simulation vector
nst = length(x0); % number of states
ninp = length(u0); % number of inputs
ntime = size(t,2); % number of time samples
if isempty(varargin)
    flag = 0;
else
    flag = varargin{1};
end
u = zeros(ninp,ntime);
[row_start,col_start] = find(abs(t-tu(:,1)) < dt^3);
[row_end,col_end] = find(abs(t-tu(:,2)) < dt^3);
id_ustart = col_start(row_start);
id_uend = col_end(row_end);
u(:,col_start(row_start):col_end(row_end)) = repmat(u_per,1,length(t...
    (id_ustart:id_uend)));
y_l = zeros(ntime, nst);
y_lq = zeros(ntime, length(ids_lq));
y_lqc = zeros(ntime, length(ids_lqc));
y_nl = zeros(ntime, nst);
tot = ntime;
%-------------------------- End of Initialization -------------------------

%------------------------ Linear System Simulation ------------------------
%cur = 0;
if flag == 1 || flag == 0
    nbytes = fprintf('Linear System Simulation...\n');
    % Initial condition:
    tmp_x = x_per;
    % Simulation procedure:
    for i = 1:ntime
        %cur = cur + 1;
        y_l(i, :) = tmp_x;
        dx = A_l * tmp_x + B_l * u(:,i);
        tmp_x = tmp_x + dx * dt;
        %step = cur / tot;
        %fprintf(repmat('\b',1,nbytes))
        %nbytes = fprintf('Linear System Simulation: %d%% \n',round(step*100));
    end
    fprintf(repmat('\b',1,nbytes))
end
%--------------------- End of Linear System Simulation --------------------

%------------------ Linear-Quadratic System Simulation --------------------
if flag == 2 || flag == 0
    nbytes = fprintf('Linear-Quadratic System Simulation...\n');
    % Initial condition:
    tmp_x = [x_per; kron(x_per,x_per)];
    tmp_x = tmp_x(ids_lq);
    % Simulation procedure:
    if ninp == 1
        for i = 1:ntime
            y_lq(i, :) = tmp_x;
            dx = A_lq * tmp_x + N_lq.N1 * tmp_x * u(1,i) + B_lq * u(:,i);
            tmp_x = tmp_x + dx * dt;
        end
    end
    
    if ninp == 2
        for i = 1:ntime
            y_lq(i, :) = tmp_x;
             dx = A_lq * tmp_x + N_lq.N1 * tmp_x * u(1,i) + N_lq.N2 * tmp_x * ...
                 u(2,i) + B_lq * u(:,i);
            tmp_x = tmp_x + dx * dt;
        end
    end
fprintf(repmat('\b',1,nbytes)) 
%     tic
%     for i = 1:ntime
%         y_lq(i, :) = tmp_x;
%         BilinTerm = 0;
%         for n = 1:ninp
%             BilinTerm = BilinTerm + eval(['N.N', num2str(n),'* tmp_x * u(n,i)']);
%         end
%         dx = A_lq * tmp_x + B_lq * u(:,i) + BilinTerm;
%         tmp_x = tmp_x + dx * dt;
%     end
%     toc
end
%--------------- End of Linear-Quadratic System Simulation ----------------

%--------------- Linear-Quadratic-Cubic System Simulation -----------------
if flag == 3 || flag == 0
    nbytes = fprintf('Linear-Quadratic-Cubic System Simulation...\n');
    % Initial condition:
    tmp_x = [x_per; kron(x_per,x_per); kron(kron(x_per,x_per),x_per)];
    tmp_x = tmp_x(ids_lqc);
    % Simulation procedure:
    if ninp == 1
        for i = 1:ntime
            y_lqc(i, :) = tmp_x;
            dx = A_lqc * tmp_x + N_lqc.N1 * tmp_x * u(1,i) + B_lqc * u(:,i);
            tmp_x = tmp_x + dx * dt;
        end
    end
    
    if ninp == 2
        for i = 1:ntime
            y_lqc(i, :) = tmp_x;
             dx = A_lqc * tmp_x + N_lqc.N1 * tmp_x * u(1,i) + N_lqc.N2 ...
                 * tmp_x * u(2,i) + B_lqc * u(:,i);
            tmp_x = tmp_x + dx * dt;
        end
    end
fprintf(repmat('\b',1,nbytes))
end
%------------ End of Linear-Quadratic-Cubic System Simulation -------------

%--------------------- Nonlinear System Simulation ------------------------
if flag == 4 || flag == 0
    nbytes = fprintf('Nonlinear System Simulation...\n');
    % Initial conditions:
    u_nl = u0 + u;
    tmp_x = x0 + x_per;
    % Simulation procedure:
    for i = 1:ntime
        y_nl(i, :) = tmp_x;
        dx = fsolvseq(sys_eq,tmp_x,u_nl(:,i));
        tmp_x = tmp_x + dx * dt;
    end
    % %opt = odeset('MaxStep',dt);
    % tic
    % [t_nl,y_nl] = ode15s(@(t,x) sys_eq(t,x,u_nl),tspan,tmp_x);
    % toc
    fprintf(repmat('\b',1,nbytes))
end
%------------------ End of Nonlinear System Simulation --------------------

%-------------------------------- Plotting --------------------------------
% if flag == 1 || flag == 0
%     figure()
%     plot(t,y_l)
%     title('linear-states')
% end
% if flag == 2 || flag == 0
%     figure()
%     %plot(t,y_lq(:,1:nst))
%     plot(t,y_lq)
%     title('linear-quadratic-states')
% end
% if flag == 3 || flag == 0
%     figure()
%     plot(t,y_nl)
%     title('nonlinear-states (absolute values)')
%     figure()
%     plot(t,y_nl-x0.')
%     title('nonlinear-states (deviation from x_0)')
% end
% if flag == 0
%     figure()
%     plot(t,u)
%     title('inputs')
%     figure()
%     plot(t,u_nl)
%     title('inputs-nonlin')
% end
%----------------------------- End of Plotting ----------------------------
end

%=============================== FUNCTIONS ================================

%-----------------------SOLVE the DYNAMIC EQUATIONS------------------------
function [fx] = fsolvseq(eq,x,u)
    fx = zeros(size(eq));
    for k = 1:size(eq,1)
        fx(k)  = eq{k,1}(x,u);   
    end
end
%--------------------END of SOLVE the DINAMIC EQUATIONS--------------------