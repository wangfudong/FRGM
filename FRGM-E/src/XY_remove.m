function [Y,rest,removerate, rest_out] = XY_remove(X,Y,order,remove_iter,type)

if nargin < 4 || isempty(remove_iter)
    remove_iter = 10;
end
if nargin < 5 || isempty(type)
    type = repmat([1,0],1,remove_iter/2);
end

Y_tmp = Y;
LX = length(X(:,1));
LY = length(Y_tmp(:,1));
rest_out = zeros(remove_iter,LY);
removerate = zeros(1,remove_iter);
removerate(1) = LX/LY;
rest_out(1,:) = 1:LY;

func_non = '1.21';
func_con = '1.21'; 
option.lam_convex = 1;
option.lam_nonvex = 0.1;


SX1 = M_points(X,X);
SX1 = 1./(SX1 + 1*(SX1==0)).*(SX1 > 0);
SX1 = SX1/max(SX1(:));
SX2 = neighbor1(X);
SX2 = SX2 + SX2';
SX2 = SX1.*SX2;

%SX2 = SX1;

option.tol_cg = 1.0e-7;
option.Mexist = 0.1;
option.maxiter_nonvex = 100;
option.maxiter_convex = 100;
option.unary = 1;
option.lam_sparse = 0;

cnt = 1;

while cnt < remove_iter
    
    M = M_shape(X,Y,1/8,2,0);
    
    Map_ini = asgHun(-M);
    
    option.convex = type(cnt);
    
    switch type(cnt)
        case 1
            option.geofunc = func_con;
            if cnt > 1
                SX = SX2;
            else
                SX = SX1;
            end
        case 0
            option.geofunc = func_non;
            SX = SX1;
    end
    
    [Map] = FW_ROT(X,Y,M,SX,Map_ini,option);
    
    [~,rest] = Y_removed(X,Y,Y_tmp,Map);
    Y = Y_tmp(rest,:);
    cnt = cnt + 1;
    removerate(cnt) = rateremove(order,rest);
    rest_out(cnt,1:length(rest)) = rest;
end

