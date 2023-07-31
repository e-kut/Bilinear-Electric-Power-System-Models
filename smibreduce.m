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
        apord = p; % approximation order
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
rr = 0;
rrr = 0;
ID = {};
exception = [];
for p = 1:size(A,1)
    if find(exception == p) ~= 0
        continue
    end
    s = s + 1;
    id(s,1) = p; % list of variable numbers without repetition
    tmp = find(X(p) == X);
    ID{p,1} = tmp;
    ltmp = length(tmp);
    if ltmp > 1
        exception = [exception; tmp(2:end)];
        if ltmp == 2
          r = r + 1;
          rmv(r,1) = tmp(2); % list of twice duplicate variable numbers to be deleted
          svd(r,1) = tmp(1); % list of twice duplicate variable numbers to be saved
        elseif ltmp == 3
            rr = rr + 1;
            rmv2_2(rr,1) = tmp(2); % list of three times duplicate variable numbers to be deleted
            rmv2_3(rr,1) = tmp(3); % list of three times duplicate variable numbers to be deleted
            svd2(rr,1) = tmp(1); % list of three times duplicate variable numbers to be saved
        elseif ltmp == 6
            rrr = rrr + 1;
            rmv3_2(rrr,1) = tmp(2);
            rmv3_3(rrr,1) = tmp(3); % list of three times duplicate variable numbers to be deleted
            rmv3_4(rrr,1) = tmp(4); % list of three times duplicate variable numbers to be deleted
            rmv3_5(rrr,1) = tmp(5); % list of three times duplicate variable numbers to be deleted
            rmv3_6(rrr,1) = tmp(6); % list of three times duplicate variable numbers to be deleted
            svd3(rrr,1) = tmp(1); % list of three times duplicate variable numbers to be saved
        end
    end
end
% Reducing of the matrices
A_red = A;
A_red(:,svd) = A(:,svd) + A(:,rmv);
if rr ~=0
    A_red(:,svd2) = A(:,svd2) + A(:,rmv2_2) + A(:,rmv2_3);
    A_red(:,svd3) = A(:,svd3) + A(:,rmv3_2) + A(:,rmv3_3) + + A(:,rmv3_4) ...
        + A(:,rmv3_5) + A(:,rmv3_6);
end
A_red = A_red(id,:);
A_red = A_red(:,id);
B_red = B(id,:);
X_red = X(id);
N_red = struct;
for k = 1:size(B,2)
    name_field = string(['N' num2str(k)]);
    tN = getfield(N,name_field);
    tN(:,svd) = tN(:,svd) + tN(:,rmv);
    if rr ~= 0
        tN(:,svd2) = tN(:,svd2) + tN(:,rmv2_2) + tN(:,rmv2_3);
        tN(:,svd3) = tN(:,svd3) + tN(:,rmv3_2) + tN(:,rmv3_3) + tN(:,rmv3_4)...
            + tN(:,rmv3_5) + tN(:,rmv3_6);
    end
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