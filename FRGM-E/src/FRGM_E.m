function asgATGM = FRGM_E(X,Y,order,option)

LX = length(X(:,1));
LY = length(Y(:,1));

%------------compute the adjacency matrix SX---------%
SX1 = sqrt(bsxfun(@minus,X(:,1),X(:,1)').^2 + bsxfun(@minus,X(:,2),X(:,2)').^2);
SX1 = 1./(SX1 + 1*(SX1==0)).*(SX1 > 0);
SX1 = SX1/max(SX1(:));

if option.full == 0
    [ntriX] = delaunay(X(:,1),X(:,2));
    SX11 = zeros(LX,LX);
    for i = 1:length(ntriX(:,1))
        tri3 = sort(ntriX(i,:),'ascend');
        SX11(tri3(1),tri3(2)) = 1;
        SX11(tri3(1),tri3(3)) = 1;
        SX11(tri3(2),tri3(3)) = 1;
    end
    SX11 = SX11 + SX11';
    SX11 = (SX11>0);
    SX1 = SX1.*SX11;
end

SX2 = neighbor_del(X);
SX2 = SX2 + SX2';
SX2 = (SX2>0);
SX2 = SX1.*SX2;
  
%-----------calculate the correspondence between graphs X and Y---------------%
rtime1 = clock;
Map = FRGM_nore(X,Y,SX1,SX2,option);
rest = 1:LY;
rtime2 = clock;

asgATGM.tim = etime(rtime2,rtime1);
asgATGM.X_con = Map;
[Maprate,~] = ratemap(Map,order,rest,0);
asgATGM.X = asgHun(Map);% post-discretization
asgATGM.acc = Maprate;
