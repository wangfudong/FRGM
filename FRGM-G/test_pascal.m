%% experiments on pascal dataset

clear variables;
prSet(1);
%%
remove = 1;
removal_time = 10;
dataset = ['Carss';'Motor'];
dataL = [30,20];

datasets = 4;%3=Carss,4=Motors
inds = 18;%index of graph pairs

[Pts,Ims,nF] = fileload(datasets,inds);
X = Pts{1,1};
Y = Pts{1,2};
I1 = Ims{1,1};
I2 = Ims{1,2};

outliers = 20;
LX = nF;
LY_all = length(Y(:,1));
LY = min(LX + outliers,LY_all);
order = randperm(LX);
%order = 1:LX;

max_size1 = max(size(rgb2gray(I1)));
max_size2 = max(size(rgb2gray(I2)));
max_size = max(max_size1,max_size2);

XX = X(order,:)/max_size;
if outliers > 0
    F2_rest = Y(LX+1:end,:);
    out_inds = randperm(LY_all-LX);
    out_inds = out_inds(1:LY-LX);
    Y_out = F2_rest(out_inds,:);% add outliers ramdonly
    Y_out = F2_rest(1:LY-LX,:);
    YY = [Y(1:LX,:);Y_out]/max_size;
else
    YY = Y(1:LY,:)/max_size;
end

if outliers > 0 && remove > 0
    [YY,Y_rest,removerate, ~] = XY_remove(XX,YY,order,removal_time,[]);
    LY = length(Y_rest);
    Pts{1,1} = XX'*max_size;
    Pts{1,2} = YY'*max_size;
else
    Y_rest = 1:LY;
    Pts{1,1} = XX'*max_size;
    Pts{1,2} = YY'*max_size;
end
removerate 
%% algorithm parameter
gt = zeros(LX,LY);
for i = 1:LX
    gt(i,order(i)) = 1;
end
asgT.alg = 'truth';
asgT.map = gt;asgT.X = gt;

parKnl = st('alg', 'pas_sc'); % type of affinity: only edge distance
[pars, algs] = gmPar(2);

parG = st('link', 'full'); % Delaunay triangulation for computing the graphs
gphs = newGphUs(Pts, parG);
%% Affinity
[MM,~,~] = M_shape(Pts{1,1}',Pts{1,2}',1/8,2);
gphs{1}.sc = MM*1;
[KP, KQ] = conKnlGphPQU(gphs, parKnl);%%parKnl 可以用 pas 表示Pascal 数据集
K = conKnlGphKU(KP, KQ, gphs);
Ct = ones(size(KP));

%% undirected graph -> directed graph (for FGM-D)
gphDs = gphU2Ds(gphs);
KQD = [KQ, KQ; KQ, KQ];

%% GA
asgGa = gm(K, Ct, asgT, pars{1}{:});
asgGa.acc = ratemap(asgHun(asgGa.X),order,Y_rest,0);

%% PM
asgPm = pm(K, KQ, gphs, asgT);
asgPm.acc = ratemap(asgHun(asgPm.X),order,Y_rest,0);

%% SM
asgSm = gm(K, Ct, asgT, pars{3}{:});
asgSm.acc = ratemap(asgHun(asgSm.X),order,Y_rest,0);

%% SMAC
asgSmac = gm(K, Ct, asgT, pars{4}{:});
asgSmac.acc = ratemap(asgHun(asgSmac.X),order,Y_rest,0);

%% IPFP-U
asgIpfpU = gm(K, Ct, asgT, pars{5}{:});
asgIpfpU.acc = ratemap(asgHun(asgIpfpU.X),order,Y_rest,0);

%% IPFP-S
asgIpfpS = gm(K, Ct, asgT, pars{6}{:});
asgIpfpS.acc = ratemap(asgHun(asgIpfpS.X),order,Y_rest,0);

%% RRWM
asgRrwm = gm(K, Ct, asgT, pars{7}{:});
asgRrwm.acc = ratemap(asgHun(asgRrwm.X),order,Y_rest,0);

%% FGM-U
asgFgmU = fgmU(KP, KQ, Ct, gphs, asgT, pars{8}{:});
asgFgmU.acc = ratemap(asgHun(asgFgmU.X),order,Y_rest,0);

%% FGM-D
asgFgmD = fgmD(KP, KQD, Ct, gphDs, asgT, pars{9}{:});
asgFgmD.acc = ratemap(asgHun(asgFgmD.X),order,Y_rest,0);

%% FRGM
option.M_exist = 0.01;
option.alpha2 = 1;
option.maxiter =100;
option.active = 1;
option.q_norm = 1;
option.type = 'inner';

M = M_shape(XX,YY,1/8,2,0);

DX = M_points(XX,XX);
DY = M_points(YY,YY);

DXX = DX/max(DX(:));
DYY = DY/max(DY(:));
SX = 1./(DXX + 1*(DXX==0)).*(DXX > 0);
SX = SX/max(SX(:));
%SX = (SX>0);

%Map_ini = asgHun(-M);
Map_ini = sinkhorn_OT(ones(LX,1),ones(LY,1)*(LX/LY),M,1/200,1.0e-7,5000,0);
sigx = std(DXX(DXX(:)>0))^2;
sigy = std(DYY(DYY(:)>0))^2;
sig = [sigx,sigy];
weight = [1,1];


asgFRGM = FRGM_Gen(Map_ini,M,DXX,DYY,SX,option,asgT,sig,weight);
asgFRGM.rate1 = ratemap(asgHun(asgFRGM.map1),order,Y_rest,0);
asgFRGM.rate2 = ratemap(asgHun(asgFRGM.map2),order,Y_rest,0);
%% MPM
asgMPM = MPM_GM(X,Y(Y_rest,:),order,K,2500);
asgMPM.acc = ratemap(asgHun(asgMPM.X),order,Y_rest,0);
