function [asgFRGM] = FRGM_testing(Map_ini,M,WX,WY,SX,option,asgt)
switch option.type
    case 'wa'
        [map,kpout] = GM_Wasserstein(Map_ini,M,WX,WY,SX,option);
    case 'inner'
        if ~isempty(asgt.sig_factor)
                        sigx = asgt.sig_factor*std(WX(WX(:)>0)).^2;
                        sigy = asgt.sig_factor*std(WY(WY(:)>0)).^2;

        else
                        sigx = 1*std(WX(WX(:)>0)).^2;
                        sigy = 1*std(WY(WY(:)>0)).^2;
            
            %             sigx = 1*median(WX(WX(:)>0)).^2;
            %             sigy = 1*median(WY(WY(:)>0)).^2;

        end
        %sigy = sigx;
        %sigx = 0.5^2;sigy = sigx;
        
        WXX = exp(-WX.^2/sigx);WYY = exp(-WY.^2/sigy);
%         [ux,vx] = eig(exp(-WX.^2/sigx));
%         [uy,vy] = eig(exp(-WX.^2/sigy));
%         WXX = exp(-WX.^2/sigx)/vx;
%         WYY = exp(-WY.^2/sigy)/vy;
        [map,kpout] = GM_PathorDis(Map_ini,M,WXX,WYY,SX,option);
    case 'metric'
        option.func = 'inner';
        [map,kpout] = GM_PathorDis(Map_ini,M,WX,WY,SX,option);
    case 'dis'
        [map,kpout] = GM_PathorDis(Map_ini,M,WX,WYSX,option);
end

LX = length(M(:,1));
LY = length(M(1,:));
right = 0;
asg = asgHun(map);
for i = 1:LX
    ind1 = find(asg(i,:) ==1 );
    ind2 = find(asgt.map(i,:)==1);
    right = right + (ind1==ind2);
end
asgFRGM.rate1 = right/LX;
asgFRGM.kpout1 = kpout;
asgFRGM.map1 = map;

%base = eye(LY);

% switch option.type
%     case 'inner'
%         M2 = zeros(LX,LY);
%         for i = 1:LX
%             for j = 1:LY
%                 M2(i,j) = (map(i,:)*WYY*map(i,:)'+1-2*(map(i,:)*WYY(:,j))).^(1);
%             end
%         end
%         %Dis_Y = (2*map*(1-WYY)).^(0.5);
%     case 'metric'
%         M2 = map*WY;
% end

%D = inner_metric(WYY,map,eye(LY));

M2 = map*WY;
%M2 = map*Dis_Y;

%M22 = map*WY;

%Map_ini = map;
option.M_exist = 1;
option.alpha2 = 1.0e0;
Map_ini = asgHun(-M2);
%[Map_ini] = sinkhorn_OT(ones(LX,1),ones(LY,1)*(LX/LY),M,1/200,1.0e-7,5000,0);
switch option.type
    case 'wa'
        [map_wa,kpout] = GM_Wasserstein(Map_ini,M2,WX,WY,SX,option);
    case 'inner'
        WX_maped = map*WYY*map';
        [map_wa,kpout] = GM_PathorDis(Map_ini,M2,WX_maped,WYY,SX,option);
    case 'metric'
        option.func = 'inner';
        WX_maped = map*WY*map';
        [map_wa,kpout] = GM_PathorDis(Map_ini,M2,WX_maped,WY,SX,option);
    case 'dis'
        [map_wa,kpout] = GM_PathorDis(Map_ini,M2,WX,WY,SX,option);
end

%map_wa = asgHun(-M2);

right = 0;
asg = asgHun(map_wa);
for i = 1:LX
    ind1 = find(asg(i,:) ==1 );
    ind2 = find(asgt.map(i,:)==1);
    right = right + (ind1==ind2);
end

asgFRGM.rate2 = right/LX;
asgFRGM.kpout2 = kpout;
asgFRGM.map2 = map_wa;

right = 0;
asg = asgHun(-M2);
for i = 1:LX
    ind1 = find(asg(i,:) ==1 );
    ind2 = find(asgt.map(i,:)==1);
    right = right + (ind1==ind2);
end
asgFRGM.rate3 = right/LX;



