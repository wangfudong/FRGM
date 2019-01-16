function [asgFRGM] = FRGM_Gen(Map_ini,M,WX,WY,SX,option,asgt,sig,weight)
% main function for FGRM-G
% Input:
%     Map_ini:initialization
%     M:unary term, dissimilarity between two graphs
%     WX,WY:pairwise edge attributes of two graphs
%     SX: weighted adjacency matrix of graph X
%     option: parameters
%     asgt: ground-truth
%     sig:sigma of X and Y
%     weight: weights of the unary and pairwise term for the second objective function
%Output:
%       asgFRGM: accuracy and correspondence map 
switch option.type
    case 'inner'
        WXX = exp(-WX.^2/sig(1));WYY = exp(-WY.^2/sig(2));
        WXX = WXX.*((WX>0)+eye(size(WX)));WYY = WYY.*((WY>0)+eye(size(WY)));
        [map,~] = GM_PathorDis(Map_ini,M,WXX,WYY,SX,option);
    case 'metric'
        [map,~] = GM_PathorDis(Map_ini,M,WX,WY,SX,option);
end

LX = length(M(:,1));
right = 0;
asg = asgHun(map);
for i = 1:LX
    ind1 = find(asg(i,:) ==1 );
    ind2 = find(asgt.map(i,:)==1);
    right = right + (ind1==ind2);
end
asgFRGM.rate1 = right/LX;
asgFRGM.map1 = map;

max_wy = max(WY(:))*2;
M2 = map*(WY+max_wy*(WY==0)-max_wy*eye(size(WY)));
option.M_exist = weight(1);
option.alpha2 = weight(2);
Map_ini = asgHun(-M2);
%Map_ini = map;

switch option.type
    case 'inner'
        WX_maped = map*WYY*map';
        [map_wa,~] = GM_PathorDis(Map_ini,M2,WX_maped,WYY,SX,option);
    case 'metric'
        WX_maped = map*WY*map';
        [map_wa,~] = GM_PathorDis(Map_ini,M2,WX_maped,WY,SX,option);
end

right = 0;
asg = asgHun(map_wa);
for i = 1:LX
    ind1 = find(asg(i,:) ==1 );
    ind2 = find(asgt.map(i,:)==1);
    right = right + (ind1==ind2);
end
asgFRGM.rate2 = right/LX;
asgFRGM.map2 = map_wa;