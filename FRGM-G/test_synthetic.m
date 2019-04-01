%% test on synthetic data
inlier = 40;
outlier = 40;
sig = 0.03;

X = randn(inlier,2);
LX = inlier;
LY = inlier + outlier;
Yout = randn(outlier,2);
Y = [X+ sig*randn(inlier,2);Yout];

max_size = max(max(abs(X(:))),max(abs(Y(:))));

order = 1:LX;
%order = randperm(LX);
X = X(order,:);

XX = X/max_size;
YY = Y/max_size;
figure,plot(XX(:,1),XX(:,2),'r.',YY(:,1),YY(:,2),'b+')

gt = zeros(LX,LY);
for i = 1:LX
    gt(i,order(i)) = 1;
end
asgT.map = gt;

%% FRGM
option.M_exist = 0.01;%0.5
option.alpha2 = 1-0.01;%0.5
option.maxiter = 100;
option.active = 1;
option.q_norm = 1;
option.type = 'inner';

M = M_shape(XX,YY,1/8,2,0);

DX = M_points(XX,XX);
DY = M_points(YY,YY);
max_fact = max([max(DX(:)),max(DY(:))]);
DXX = DX/max_fact;
DYY = DY/max_fact;

SX = 1./(DXX + 1*(DXX==0)).*(DXX > 0);
SX = SX/max(SX(:));
%SX = (SX>0);

Map_ini = asgHun(-M);
sig = [0.5,0.5].^2;
weight = [1,1];
asgFRGM = FRGM_Gen(Map_ini,M,DXX,DYY,SX,option,asgT,sig,weight)

%% plot
map1 = asgHun(asgFRGM.map2);
%map1 = asgHun(Map_ini);
YP = YY;
YP(:,1) = YP(:,1) + 3;
YYP = map1*YP;

matched = zeros(1,length(order));
for ii = 1:length(order)
    matched1(ii) = find(map1(ii,:)==1);
    matched(ii) = (matched1(ii)==order(ii));
end
rate_end = sum(matched)/length(order)

SX1 = neighbor1(XX);SX1 = (SX1>0);
SX2 = neighbor1(YY);SX2 = (SX2>0);



figure('color','white');
plot(XX(:,1),XX(:,2),'r.','markersize',30);hold on;
for i=1:length(order)
    rep = find(SX1(i,:)==1);
    line([repmat(XX(i,1),length(rep)),XX(rep,1)]',[repmat(XX(i,2),length(rep)),XX(rep,2)]','color','k');
    hold on;
end
plot(YP(:,1),YP(:,2),'b.','markersize',20);hold on;
plot(YP(order,1),YP(order,2),'r.','markersize',30);hold on;
for i=1:(length(YP(:,1)))
    rep = find(SX2(i,:)==1);
    line([repmat(YP(i,1),length(rep)),YP(rep,1)]',[repmat(YP(i,2),length(rep)),YP(rep,2)]','color','k');
    hold on;
end
axis equal;
axis off;
line([XX(:,1),YYP(:,1)]',[XX(:,2),YYP(:,2)]','color','g','linewidth',2);
for ii = 1:length(order)
    if matched(ii)==0
        line([XX(ii,1),YYP(ii,1)],[XX(ii,2),YYP(ii,2)],'color','m','linewidth',2);hold on;
    end
end
