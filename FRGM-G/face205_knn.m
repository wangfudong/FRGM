%% test fgmdemos
%%
dataset = ['face205'];
dataL = 37;

knn_rate = 0.9;

LX = 40;
LY = 50;

knei_x = LX*knn_rate;
knei_y = LY*knn_rate;
inds = 15;
orders = [336,336+30*inds];
[DX,DY,hks1,hks2] = get_graph_face205(orders);

order = randperm(LY);
order = order(1:LX);

%% algorithm parameter
gt = zeros(LX,LY);
for i = 1:LX
    gt(i,order(i)) = 1;
end
asgT.map = gt;asgT.X = gt;

%%    FRGM
option.M_exist = 0.01;
option.alpha2 = 1;
option.maxiter = 100;
option.active = 1;
option.q_norm = 1;
option.type = 'inner';

hks11 = hks1(order,:);
M = measure_hks(hks11,hks2,'E');

DXX = DX/max(DX(:));
DYY = DY/max(DY(:));
DXX = DXX(order,:);
DXX = DXX(:,order);
SX = 1./(DXX + 1*(DXX==0)).*(DXX > 0);
SX = SX/max(SX(:));
%SX = (SX>0);

sx_knn = neighbor_knn(DXX,knei_x);
sy_knn = neighbor_knn(DYY,knei_y);
SX = sx_knn.*SX;

DXX_knn = DXX.*sx_knn;
DYY_knn = DYY.*sy_knn;


Map_ini = asgHun(-M);

sig = [0.5,0.5].^2;
weight = [1,1];
tic;
[asgFRGM] = FRGM_Gen(Map_ini,M,DXX_knn,DYY_knn,SX,option,asgT,sig,weight);
toc;
