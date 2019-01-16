function [Map_next,kpout] = GM_PathorDis(Map_ini,M0,DX,DY,SX,option)
% main function for computing the correspondence
maxiter = option.maxiter;
min_iter = 2;

ST = 20;
b = 1/2;
kp = 0;
right = 1;
cnt = 1;
map_tmp = Map_ini;

option.gra = 1;
[~,M] = func_handle(M0,Map_ini,SX,DX,DY,option);

kpout = zeros(1,option.maxiter);
LX = length(DX(:,1));
LY = length(DY(:,1));

active = option.active;
hist_direc = zeros(LX,LY,active);
weights = ones(1,active)/active;
%weights = exp(-((1:active)-active).^2/active);weights = weights/sum(weights);

condition1 = 1;
condition2 = 1;

while condition1 > 0 || condition2 > 0
    
    %      Map = asgHun(-M);% Hungarian algorithm for LAP
    
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
    
     if cnt <= active
        dmap = Map - map_tmp;% dmap=v - u
        hist_direc(:,:,cnt) = Map;
    else
        hist_direc(:,:,active+1) = Map;
        hist_direc(:,:,1) = [];
        dmap = 0;
        for ii = 1:active
            dmap = dmap + weights(ii)*(hist_direc(:,:,ii));
        end
        dmap = dmap - map_tmp;
    end
    
    left = -sum(sum(M.*dmap));
    right = 1;
    breakyn = 1;
    if abs(left) >= 1.0e-7
        while kp <= ST && abs(right) >= 1.0e-7 && breakyn > 0
            option.gra = 0;
            F1 = func_handle(M0,map_tmp,SX,DX,DY,option);
            F2 = func_handle(M0,map_tmp+b^kp*dmap,SX,DX,DY,option);
            right = F1 - F2;
            if right >= -1.0e-7
                %break;
                option.gra = 0;
                right0 = F1 - func_handle(M0,map_tmp+b^(kp+1)*dmap,SX,DX,DY,option);
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
    
    kpout(cnt) = kp;
    cnt = cnt +1;
    condition1 = (cnt <= min_iter);
    condition2 = (cnt <= maxiter)*(kp <= ST)*(abs(right) >= 1.0e-7);
    
    if  kp <= ST
        option.gra = 1;
        [~,M] = func_handle(M0,Map_next,SX,DX,DY,option);
    end
    
end

disp(['  Last cg_Iter:',num2str(cnt-1),' Last cg_Err: ',num2str(right) ' k= ' num2str(kp)]);

end





