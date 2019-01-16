function [Map_next] = FW_Min(X,Y,M,SX,Map_ini,option,approximate)
% Frank-Wolfe method-based minimization:
% min_P a1*<P,M> + a2*L(P) + a3*S(P)
%Input:
% X, Y: nodes of two graphs
% M: dissimilarity matrix between X and Y: M(Xi,Yj)
% SX: similarity or weight matrix of graph X
% Map_ini: inlitialization of P
% option: options of parameters
% approximate: whether to use AFW or FW
%Output:
% Map_next: the final correspondence map between two
%% initialization and set parameters
Mexist = option.Mexist;

if option.convex == 0
    if option.unary == 1
        maxiter = option.maxiter_nonvex;
        a2 = option.lam_nonvex;
        [Fgeo,grad_geo] = geofunc_handle(X,Y,SX,option);
        F = @(P)( sum(sum(Mexist*M.*P)) +  a2*Fgeo(P) );
        gFt = @(P)(Mexist*M + a2*grad_geo(P));
    else
        maxiter = option.maxiter_nonvex;
        a2 = option.lam_nonvex;
        [Fgeo,grad_geo] = geofunc_handle(X,Y,SX,option);
        F = @(P)( Mexist*sum(sum((X - P*Y).^2)) +  a2*Fgeo(P) );
        gFt = @(P)(Mexist*2*(P*Y-X)*Y' + a2*grad_geo(P));
        
    end
else
    if option.unary == 1
        maxiter = option.maxiter_convex;
        a2 = option.lam_convex;
        a3 = option.lam_sparse;
        [Fgeo,grad_geo] = geofunc_handle(X,Y,SX,option);
        F = @(P)( sum(sum(Mexist*M.*P)) +  a2*Fgeo(P) + a3*sum(sum(abs(P))));
        gFt = @(P)(Mexist*M + a2*grad_geo(P) + a3);
    else
        maxiter = option.maxiter_convex;
        a2 = option.lam_convex;
        a3 = option.lam_sparse;
        [Fgeo,grad_geo] = geofunc_handle(X,Y,SX,option);
        F = @(P)( Mexist*sum(sum((X - P*Y).^2)) +  a2*Fgeo(P) + a3*sum(sum(abs(P))));
        gFt = @(P)(Mexist*2*(P*Y-X)*Y' + a2*grad_geo(P) + a3);
    end
end


%%  initialization
if isempty(Map_ini) > 0
    error('No initial point!');
end

if approximate > 0
    ST = 20;
else
    ST = 15;
end

b = 1/2;
kp = 0;
right = 1;
cnt = 1;
map_tmp = Map_ini;
M = gFt(Map_ini);

LX = length(X(:,1));
LY = length(Y(:,1));

lambda = 1/500;
while cnt <= maxiter  && kp <= ST && abs(right) >= 1.0e-7
    
    if approximate > 0
        M = M-min(M(:));
        M = M/max(M(:));
        Map = sinkhorn_OT(ones(LX,1),ones(LY,1)*(LX/LY),M,lambda,1.0e-6,1000,0);
        %Map = partial_OT(ones(LX,1),ones(LY,1),M,lambda,1,1.0e-7,3000);
    else
        M1 = [M,ones(LX,10)*max(max(abs(M)))+100];
        % M1 = M;
        [rowsol] = lapjv(M1,1.0e0);
        Map = zeros(size(M));
        for i = 1:LX
            Map(i,rowsol(i)) = 1;
        end
    end
    
    % for less ierations
    if cnt <= 30
        kp = 0;
    else
        kp = 5;
    end
    
    %-----------inexact line-search---------------%
    dmap = Map - map_tmp;% dmap=v - u
    left = -sum(sum(gFt(map_tmp).*dmap));
    right = 1;
    breakyn = 1;
    if abs(left) >= 1.0e-7
        F_tmp = F(map_tmp);
        while kp <= ST && abs(right) >= 1.0e-7 && breakyn > 0
            right = F_tmp - F(map_tmp+b^kp*dmap);
            if right >= -1.0e-7
                %break;
                right0 = F_tmp - F(map_tmp+b^(kp+1)*dmap);
                if right0 > right
                    kp = kp + 1;
                    break;
                else
                    break;
                end
            else
                kp = kp + 1;
            end
        end
        Map_next = map_tmp + b.^kp*dmap;
        map_tmp = Map_next;
    else
        kp = ST + 1;
        Map_next = Map;
    end
    
    cnt = cnt +1;
    
    if  kp <= ST
        M = gFt(Map_next);
    end
    
end



