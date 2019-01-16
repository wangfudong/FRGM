function [Map_next] = FW_ROT(X,Y,M,SX,Map_ini,option)
% Conditional Gradient (Frank-Wolfe) descent optimization for:
% min <P,M> + a1*L(P) + a2*S(P)
% inputs:
% X, Y: vertex (or samplings) of two graphs
% M: ground distance matrix between X and Y: M(Xi,Yj)
% SX,SY: similarity or weight matrix of graph X and graph Y
% option: struct to store parameters

%% initialization and set parameters
Mexist = option.Mexist;

M_EOT = M;

if option.convex == 0
    if option.unary == 1
        maxiter = option.maxiter_nonvex;
        a1 = option.lam_nonvex;
        [Fgeo,grad_geo] = geofunc_handle(X,Y,SX,option);
        F = @(P)( sum(sum(Mexist*M_EOT.*P)) +  a1*Fgeo(P) );
        gFt = @(P)(Mexist*M_EOT + a1*grad_geo(P));
    else
        maxiter = option.maxiter_nonvex;
        a1 = option.lam_nonvex;
        [Fgeo,grad_geo] = geofunc_handle(X,Y,SX,option);
        F = @(P)( Mexist*sum(sum((X - P*Y).^2)) +  a1*Fgeo(P) );
        gFt = @(P)(Mexist*2*(P*Y-X)*Y' + a1*grad_geo(P));
        
    end
else
    switch option.unary
        case 1
            maxiter = option.maxiter_convex;
            a2 = option.lam_convex;
            a3 = option.lam_sparse;
            [Fgeo,grad_geo] = geofunc_handle(X,Y,SX,option);
            F = @(P)( sum(sum(Mexist*M_EOT.*P)) +  a2*Fgeo(P) + a3*sum(sum(abs(P))));
            gFt = @(P)(Mexist*M_EOT + a2*grad_geo(P) + a3);
            
        case 2
            
            maxiter = option.maxiter_convex;
            a2 = option.lam_convex;
            a3 = option.lam_sparse;
            [Fgeo,grad_geo] = geofunc_handle(X,Y,SX,option);
            F = @(P)( Mexist*sum(sum((X - P*Y).^2)) +  a2*Fgeo(P) + a3*sum(sum(abs(P))));
            gFt = @(P)(Mexist*2*(P*Y-X)*Y' + a2*grad_geo(P) + a3);
    end
    
end


%% cg optimization

if isempty(Map_ini) > 0
    error('No initial point!');
end
min_iter = 10;
ST = 15;
b = 1/2;
kp = 0;
right = 1;
cnt = 1;
map_tmp = Map_ini;
M = gFt(Map_ini);

LX = length(X(:,1));

condition1 = 1;
condition2 = 1;

while condition1 > 0 || condition2 > 0
    
%     Map = asgHun(-M);% Hungarian algorithm for LAP
    

    M1 = [M,ones(LX,0)*max(max(abs(M)))+100];
    % M1 = M;
    [rowsol] = lapjv(M1,1.0e0);% LAPJV algorithm for LAP
    Map = zeros(size(M));
    for i = 1:LX
        Map(i,rowsol(i)) = 1;
    end
    
%     [~,P] = sparseAssignmentProblemAuctionAlgorithm(sparse(max(max(M))-M));
%     Auction algorithm for LAP, much faster for equal-sized case
%     Map = full(P);

    if cnt <=30% for less iteration
        kp = 0;
    else
        kp = 5;
    end
    dmap = Map - map_tmp;% dmap=v - u
    left = -sum(sum(gFt(map_tmp).*dmap));
    right = 1;
    breakyn = 1;
    t1 = clock;
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
    condition1 = (cnt <= min_iter);
    condition2 = (cnt <= maxiter)*(kp <= ST)*(abs(right) >= 1.0e-7);
    
    if  kp <= ST
        M = gFt(Map_next);
    end

end


disp(['  Last cg_Iter:',num2str(cnt-1),' Last cg_Err: ',num2str(right) ' k= ' num2str(kp)]);

