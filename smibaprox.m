function [sys_mat, varargout] = smibaprox(sys_eq, x0, u0)
%------------------------------ Description -------------------------------
%
% This function construct linear, linear-quadratic and linear-quadratic-
% cubic systems to study the dynamics in deviations from the operating 
% point (x0,u0)
%
% INPUTS:   sys_eq  -   (cell array) set of the system equations in handle 
%                       function form
%           x0      -   (column-vector) initial steady-state values of the
%                       state variables
%           u0      -   (column-vector) initial steady-state values of the
%                       input variables
%
% OUTPUTS:  sys_mat -   structure containing A_l, B_l, A_lq, B_lq, N_lq,
%                       A_lqc, B_lqc, N_lqc matrices
%           A_l     -   dynamics matrix of the linear system
%           B_l     -   input matrix of the linear system
%           A_lq    -   dynamics matrix of the linear-quadratic system
%           B_lq    -   input matrix of the linear-quadratic system
%           N_lq    -   set of matrices of linear-quadratic bilinear terms
%           A_lqc   -   dynamics matrix of the linear-quadratic-cubic system
%           B_lqc   -   input matrix of the linear-quadratic-cubic system
%           N_lqc   -   set of matrices of linear-quadratic-cubic bilinear terms
%           ids_lq  -   numbers of saved states in the original vector of 
%                       state variables for linear-quadratic system
%           ids_lqc -   numbers of saved states in the original vector of 
%                       state variables for linear-quadratic-cubic system
%
% SYNTAX: 
%       (a) full (first, second and third order) approximation:
%           (a.a)   smibaprox(sys_eq,x,u);
%           (a.b)   [~] = smibaprox(sys_eq,x,u);
%           (a.c)   [sys_mat] = smibaprox(sys_eq,x,u);
%           (a.d)   [~,A_l,B_l,A_lq,B_lq,N_lq,A_lqc,B_lqc,N_lqc,ids_lq, ids_lqc] = 
%                   smibaprox(sys_eq,x,u);
%           (a.e)   [sys_mat,A_l,B_l,A_lq,B_lq,N_lq,A_lqc,B_lqc,N_lqc,ids_lq, ids_lqc]
%                   = smibaprox(sys_eq,x,u);
%       (b) first order approximation:
%           (b.a)   [~,A_l,B_l] = smibaprox(sys_eq,x,u);
%           (b.b)   [sys_mat,A_l,B_l] = smibaprox(sys_eq,x,u);
%       (c) first and second order approximation:
%           (c.a)   [~,A_l,B_l,A_lq,B_lq,N_lq,ids_lq] = smibaprox(sys_eq,x,u);
%           (c.b)   [sys_mat,A_l,B_l,A_lq,B_lq,N_lq,ids_lq] = smibaprox
%                   (sys_eq,x,u);
%
%
% Linear System: dx/dt = A_l * x + B_l * u
% Linear-Quadratic System: dx/dt = A_lq * x + sum_i(Ni*x*ui) + B_lq * u
%--------------------------- End of Description ---------------------------

%----------------------------- Main Procedure -----------------------------
% Function of solving overall EPS equations
F = @(x,u) fsolvseq(sys_eq,x,u);
% Jacobian
[A1x,A1u] = fjaco(F,x0,u0);
% State-space matrix of the linear system
A_l = A1x;
B_l = A1u;
% Outputs
varargout{1} = A_l;
varargout{2} = B_l;
sys_mat.A_l = A_l;
sys_mat.B_l = B_l;

if nargout ~= 3
    % Hessian
    [A2x,~,~] = fhess(F,x0,u0);
    % State-space matrix of the linear-quadratic system
    nx = size(A1x,1);
    n2x = nx + nx^2;
    nu = size(B_l,2);
    I = eye(nx, nx);
    A21 = kron(A1x, I) + kron(I, A1x);
    A_lq = zeros(n2x, n2x);
    A_lq(1:nx, 1:nx) = A1x;
    A_lq(1:nx, nx+1:n2x) = A2x;
    A_lq(nx+1:n2x, nx+1:n2x) = A21;
    B_lq = zeros(n2x, size(u0,1));
    B_lq(1:nx,:) = B_l;
    for k = 1:nu
        eval(['N_lq.N',num2str(k), '= zeros(size(A_lq));'])
        eval(['N_lq.N',num2str(k), '(nx+1:n2x, 1:nx) = kron(B_l(:,k), I)',...
            '+ kron(I, B_l(:,k));'])
    end
    
    varargout{3} = A_lq;
    varargout{4} = B_lq;
    varargout{5} = N_lq;
    varargout{6} = 1:n2x; % ids_lq
    sys_mat.A_lq = A_lq;
    sys_mat.B_lq = B_lq;
    sys_mat.N_lq = N_lq;
    sys_mat.ids_lq = varargout{6};
end

if nargout ~= 3 && nargout ~= 7
    % Cubic Hessian
    [A3x,~,~] = fchess(F,x0,u0);
    % State-space matrix of the linear-quadratic-cubic system
    n3x = n2x + nx^3;
    A22 = kron(A2x, I) + kron(I, A2x);
    A31 = kron(kron(A1x, I),I) + kron(kron(I,A1x),I) + kron(kron(I,I),A1x);
    A_lqc = zeros(n3x, n3x);
    A_lqc(1:nx, 1:nx) = A1x;
    A_lqc(1:nx, nx+1:n2x) = A2x;
    A_lqc(1:nx, n2x+1:n3x) = A3x;
    A_lqc(nx+1:n2x, nx+1:n2x) = A21;
    A_lqc(nx+1:n2x, n2x+1:n3x) = A22;
    A_lqc(n2x+1:n3x, n2x+1:n3x) = A31;
    B_lqc = zeros(n3x, size(u0,1));
    B_lqc(1:nx,:) = B_l;
    for k = 1:nu
        eval(['N_lqc.N',num2str(k), '= zeros(size(A_lqc));'])
        eval(['N_lqc.N',num2str(k), '(1:n2x, 1:n2x) = N_lq.N',num2str(k),';'])
        eval(['N_lqc.N',num2str(k), '(n2x+1:n3x, nx+1:n2x) = ',...
            'kron(kron(B_l(:,k), I), I) + kron(kron(I, B_l(:,k)), I)', ...
            '+ kron(kron(I, I), B_l(:,k));'])
    end
    varargout{6} = A_lqc;
    varargout{7} = B_lqc;
    varargout{8} = N_lqc;
    varargout{9} = 1:n2x; % ids_lq
    varargout{10} = 1:n3x; % ids_lqc
    sys_mat.A_lqc = A_lqc;
    sys_mat.B_lqc = B_lqc;
    sys_mat.N_lqc = N_lqc;
    sys_mat.ids_lq = varargout{9};
    sys_mat.ids_lqc = varargout{10};
end
%-------------------------- End of Main Procedure -------------------------
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

%-----------------------CALCULATING the JACOBIAN---------------------------
function [jacx,jacv] = fjaco(F,x,u)

    d = 1e-6;
    for k = 1:size(x,1)
        tx1 = x;
        tx2 = x;
        tx1(k) = tx1(k) + d/2;
        tx2(k) = tx2(k) - d/2;
        jacx(:,k) = (F(tx1,u) - F(tx2,u)) / d;
    end
    
    for k = 1:size(u,1)
        tv1 = u;
        tv2 = u;
        tv1(k) = tv1(k) + d/2;
        tv2(k) = tv2(k) - d/2;
        jacv(:,k) = (F(x,tv1) - F(x,tv2)) / d;
    end
end
%--------------------END of CALCULATING the JACOBIAN-----------------------

%-----------------------CALCULATING the HESSIAN----------------------------
function [hesx,hesu,hesux] = fhess(F,x,u)
    fxu = F(x,u);
    % ƒифференцирование по переменным состо€ни€ (вектор x)
    d = 1e-5;
    ast = 1:size(x,1);
    nst = size(x,1);
    for i = ast% k = 3:size(x,1) если нужно удалить первые два уравнени€
        xip = x;
        xim = x;
        xip(i) = xip(i) + d;
        xim(i) = xim(i) - d;
        fxip = F(xip,u);
        fxim = F(xim,u);
        hesx(:,(i-1)*nst+i)=(fxip - 2 * fxu + fxim) / (d^2);
        for j = ast([1:i-1 i+1:end])% зеркальные коэффициенты не удал€ютс€
            xipjp = xip;
            ximjp = xim;
            xipjm = xip;
            ximjm = xim;
            xipjp(j) = xip(j) + d;
            ximjp(j) = xim(j) + d;
            xipjm(j) = xip(j) - d;
            ximjm(j) = xim(j) - d;
            fxipjp = F(xipjp,u);
            fximjp = F(ximjp,u);
            fxipjm = F(xipjm,u);
            fximjm = F(ximjm,u);
            hesx(:,j+nst*(i-1)) = (fxipjp - fximjp - fxipjm + fximjm) /...
                (4*d^2);
        end
    end
    
    % ƒифференцирование по входам (вектор u)
    %d = 1e-2;
    ainp = 1:size(u,1);
    ninp = size(u,1);
    for i = ainp% k = 3:size(x,1) если нужно удалить первые два уравнени€
        uip = u;
        uim = u;
        uip(i) = uip(i) + d;
        uim(i) = uim(i) - d;
        fuip = F(x,uip);
        fuim = F(x,uim);
        hesu(:,(i-1)*ninp+i) = (fuip - 2 * fxu + fuim) / (d^2);
        for j = ainp([1:i-1 i+1:end])% зеркальные коэффициенты не удал€ютс€
            uipjp = uip;
            uimjp = uim;
            uipjm = uip;
            uimjm = uim;
            uipjp(j) = uipjp(j) + d;
            uimjp(j) = uimjp(j) + d;
            uipjm(j) = uipjm(j) - d;
            uimjm(j) = uimjm(j) - d;
            fuipjp = F(x,uipjp);
            fuimjp = F(x,uimjp);
            fuipjm = F(x,uipjm);
            fuimjm = F(x,uimjm);
            hesu(:,j+ninp*(i-1)) = (fuipjp - fuimjp - fuipjm + fuimjm) ...
                / (4*d^2);
        end
    end
    
    % —мешанные производные по x и u 
    %d = 1e-4;
    for i = ainp% k = 3:size(x,1) если нужно удалить первые два уравнени€
        uip = u;
        uim = u;
        uip(i) = uip(i) + d;
        uim(i) = uim(i) - d;
        for j = 1:size(x,1)
            xjp = x;
            xjm = x;
            xjp(j) = xjp(j) + d;
            xjm(j) = xjm(j) - d;
            fuipxjp = F(xjp,uip);
            fuimxjp = F(xjp,uim);
            fuipxjm = F(xjm,uip);
            fuimxjm = F(xjm,uim);
            hesux(:,j+nst*(i-1)) = (fuipxjp - fuimxjp - fuipxjm + fuimxjm)...
                / (4*d^2);
        end
    end
   
    hesx = hesx / 2;
    hesu = hesu / 2;
    hesux = hesux / 2;
    
end
%--------------------END of CALCULATING the HESSIAN------------------------

%--------------------CALCULATING the CUBIC HESSIAN-------------------------
function [hesx,hesu,hesux] = fchess(F,x,u)
    % ƒифференцирование по переменным состо€ни€ (вектор x)
    d = 1e-2;
    id = 1:size(x,1);
    nst = size(x,1);
    for i = id% k = 3:size(x,1) если нужно удалить первые два уравнени€
        idj = id;
        idj(i) = []; % набор индексов без i
        bid = (i-1)*nst^2; % начальный индекс i-го блока
        subid = (i-1)*nst; % начальный индекс i-го подблока
        xip = x;
        xip2 = x;
        xim = x;
        xim2 = x;
        xip(i) = xip(i) + d;
        xim(i) = xim(i) - d;
        xip2(i) = xip2(i) + 2*d;
        xim2(i) = xim2(i) - 2*d;
        fxip = F(xip,u);
        fxim = F(xim,u);
        fxip2 = F(xip2,u);
        fxim2 = F(xim2,u);
        hesx(:,(i-1)*(nst^2+nst)+i) = (fxip2 - 2 * fxip + 2 * fxim - ...
            fxim2) / (2 * d^3);
        for j = idj% зеркальные коэффициенты не удал€ютс€
            idk = id;
            idk(j) = []; % набор индексов без i и j
            subjd = (j-1)*nst; % начальный индекс j-го подблока
            xipjp = xip;
            ximjp = xim;
            xipjm = xip;
            ximjm = xim;
            xjp = x;
            xjm = x;
            xipjp(j) = xip(j) + d;
            ximjp(j) = xim(j) + d;
            xipjm(j) = xip(j) - d;
            ximjm(j) = xim(j) - d;
            xjp(j) = xjp(j) + d;
            xjm(j) = xjm(j) - d;
            fxipjp = F(xipjp,u);
            fximjp = F(ximjp,u);
            fxipjm = F(xipjm,u);
            fximjm = F(ximjm,u);
            fxjp = F(xjp,u);
            fxjm = F(xjm,u);
            hesx(:,bid+(i-1)*nst+j) = (fxipjp - 2*fxjp + fximjp - fxipjm ...
                + 2*fxjm - fximjm) / (2*d^3);
            hesx(:,bid+(j-1)*nst+i) = hesx(:,bid+(i-1)*nst+j);
            hesx(:,(j-1)*nst^2 +subid+i) = hesx(:,bid+(i-1)*nst+j);
            for k = idk
                xipjpkp = xipjp;
                xipjpkp(k) = xipjp(k) + d;
                ximjpkp = ximjp;
                ximjpkp(k) = ximjp(k) + d;
                xipjmkp = xipjm;
                xipjmkp(k) = xipjm(k) + d;
                ximjmkp = ximjm;
                ximjmkp(k) = ximjm(k) + d;
                xipjpkm = xipjp;
                xipjpkm(k) = xipjp(k) - d;
                ximjpkm = ximjp;
                ximjpkm(k) = ximjp(k) - d;
                xipjmkm = xipjm;
                xipjmkm(k) = xipjm(k) - d;
                ximjmkm = ximjm;
                ximjmkm(k) = ximjm(k) - d;
                fxipjpkp = F(xipjpkp,u);
                fximjpkp = F(ximjpkp,u);
                fxipjmkp = F(xipjmkp,u);
                fximjmkp = F(ximjmkp,u);
                fxipjpkm = F(xipjpkm,u);
                fximjpkm = F(ximjpkm,u);
                fxipjmkm = F(xipjmkm,u);
                fximjmkm = F(ximjmkm,u);
                hesx(:,bid+subjd+k) = (fxipjpkp - fximjpkp - fxipjmkp ...
                    + fximjmkp - fxipjpkm + fximjpkm + fxipjmkm - fximjmkm) ...
                    / (8*d^3);
            end
        end
    end
    
    % ƒифференцирование по входам (вектор u)
    d = 1e-2;
    id = 1:size(u,1);
    nst = size(u,1);
    for i = id% k = 3:size(x,1) если нужно удалить первые два уравнени€
        idj = id;
        idj(i) = []; % набор индексов без i
        bid = (i-1)*nst^2; % начальный индекс i-го блока
        subid = (i-1)*nst; % начальный индекс i-го подблока
        uip = u;
        uip2 = u;
        uim = u;
        uim2 = u;
        uip(i) = uip(i) + d;
        uim(i) = uim(i) - d;
        uip2(i) = uip2(i) + 2*d;
        uim2(i) = uim2(i) - 2*d;
        fuip = F(x,uip);
        fuim = F(x,uim);
        fuip2 = F(x,uip2);
        fuim2 = F(x,uim2);
        hesu(:,(i-1)*(nst^2+nst)+i) = (fuip2 - 2 * fuip + 2 * fuim - ...
            fuim2) / (2 * d^3);
        for j = idj% зеркальные коэффициенты не удал€ютс€
            idk = id;
            idk(j) = []; % набор индексов без i и j
            subjd = (j-1)*nst; % начальный индекс j-го подблока
            uipjp = uip;
            uimjp = uim;
            uipjm = uip;
            uimjm = uim;
            ujp = u;
            ujm = u;
            uipjp(j) = uip(j) + d;
            uimjp(j) = uim(j) + d;
            uipjm(j) = uip(j) - d;
            uimjm(j) = uim(j) - d;
            ujp(j) = ujp(j) + d;
            ujm(j) = ujm(j) - d;
            fuipjp = F(x,uipjp);
            fuimjp = F(x,uimjp);
            fuipjm = F(x,uipjm);
            fuimjm = F(x,uimjm);
            fujp = F(x,ujp);
            fujm = F(x,ujm);
            hesu(:,bid+(i-1)*nst+j) = (fuipjp - 2*fujp + fuimjp - fuipjm ...
                + 2*fujm - fuimjm) / (2*d^3);
            hesu(:,bid+(j-1)*nst+i) = hesu(:,bid+(i-1)*nst+j);
            hesu(:,(j-1)*nst^2 +subid+i) = hesu(:,bid+(i-1)*nst+j);
            for k = idk
                uipjpkp = uipjp;
                uipjpkp(k) = uipjp(k) + d;
                uimjpkp = uimjp;
                uimjpkp(k) = uimjp(k) + d;
                uipjmkp = uipjm;
                uipjmkp(k) = uipjm(k) + d;
                uimjmkp = uimjm;
                uimjmkp(k) = uimjm(k) + d;
                uipjpkm = uipjp;
                uipjpkm(k) = uipjp(k) - d;
                uimjpkm = uimjp;
                uimjpkm(k) = uimjp(k) - d;
                uipjmkm = uipjm;
                uipjmkm(k) = uipjm(k) - d;
                uimjmkm = uimjm;
                uimjmkm(k) = uimjm(k) - d;
                fuipjpkp = F(x,uipjpkp);
                fuimjpkp = F(x,uimjpkp);
                fuipjmkp = F(x,uipjmkp);
                fuimjmkp = F(x,uimjmkp);
                fuipjpkm = F(x,uipjpkm);
                fuimjpkm = F(x,uimjpkm);
                fuipjmkm = F(x,uipjmkm);
                fuimjmkm = F(x,uimjmkm);
                hesu(:,bid+subjd+k) = (fuipjpkp - fuimjpkp - fuipjmkp ...
                    + fuimjmkp - fuipjpkm + fuimjpkm + fuipjmkm - fuimjmkm) ...
                    / (8*d^3);
            end
        end
    end
%     
%     % —мешанные производные по x и u 
%     %d = 1e-4;
%     for i = ainp% k = 3:size(x,1) если нужно удалить первые два уравнени€
%         uip = u;
%         uim = u;
%         uip(i) = uip(i) + d;
%         uim(i) = uim(i) - d;
%         for j = 1:size(x,1)
%             xjp = x;
%             xjm = x;
%             xjp(j) = xjp(j) + d;
%             xjm(j) = xjm(j) - d;
%             fuipxjp = F(xjp,v,uip);
%             fuimxjp = F(xjp,v,uim);
%             fuipxjm = F(xjm,v,uip);
%             fuimxjm = F(xjm,v,uim);
%             hesux(:,j+nst*(i-1)) = (fuipxjp - fuimxjp - fuipxjm + fuimxjm)...
%                 / (4*d^2);
%         end
%     end
%     else
%         hesu = [];
%         hesux = [];
%     end


    hesx = hesx / 6;
    hesu = hesu / 6;
    hesux = [];
end
%-----------------END of CALCULATING the CUBIC HESSIAN---------------------