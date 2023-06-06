function [A_red, B_red, N_red, ids] = smibreduce (A,B,N,order)
%------------------------------ Description -------------------------------
%
% This function performs the removal of repeated states in the dynamics 
% matrix of a linear-quadratic or linear-quadratic-cubic systems
%
% INPUTS:   A       -   dynamics matrix of the unreduced system
%           B       -   input matrix of the unreduced system
%           N       -   set of matrices of bilinear terms
%           order   -   order of the linear system (initial state space 
%                       dimension)
%
% OUTPUTS:  A_red   -   dynamics matrix of the reduced system
%           B_red   -   input matrix of the reduced system
%           N_red   -   set of reduced matrices of bilinear terms
%           ids     -   numbers of saved states in the original vector of 
%                       state variables
%
%--------------------------- End of Description ---------------------------

%----------------------------- Main Procedure -----------------------------
% Determine the approximation order
tmp = 0;
for p = 1:3
    tmp = tmp + order^p;
    ind = tmp / size(A,1);
    if ind == 1
        apord = p;
        break
    end
end
% Initializing of the states vector
syms x [order,1]
tmpx = x;
X = x;
for n = 1:apord-1
    tmpx = kron(tmpx,x);
    X = [X;tmpx];
end
% Finding repetitive states
s = 0;
r = 0;
ID = {};
for p = 1:size(A,1)
    tmp = find(X(p) == X);
    ID{p,1} = tmp;
    if tmp(1) == p
        s = s + 1;
        id(s,1) = p; % list of variable numbers without repetition
    else
        r = r + 1;
        rmv(r,1) = p; % list of duplicate variable numbers to be deleted
        svd(r,1) = tmp(1); % list of duplicate variable numbers to be saved
    end
end
% Reducing of the matrices
A_red = A;
A_red(:,svd) = A(:,svd) + A(:,rmv);
A_red = A_red(id,:);
A_red = A_red(:,id);
B_red = B(id,:);
X_red = X(id);
N_red = struct;
for k = 1:size(B,2)
    name_field = string(['N' num2str(k)]);
    tN = getfield(N,name_field);
    tN(:,svd) = tN(:,svd) + tN(:,rmv);
    Nred = tN(id,:);
    Nred = Nred(:,id);
    N_red = setfield(N_red,name_field,Nred);
end
ids = id;
% reduce_sys.A = A_red;
% reduce_sys.A1 = A_red(1:24,1:24);
% reduce_sys.B = B_red;
% reduce_sys.B1 = B_red(1:24,:);
% reduce_sys.N = N;
% reduce_sys.x_full = X_red;
%-------------------------- End of Main Procedure -------------------------
end