%% test cmu house/hotel
clear variables;
prSet(1);
%%
dataset = 1;% 1 for house (111 images), 2 for hotel (101 images)

LX = 20;%LX<=30;
LY = 30;%LY<=30;

inds = [1,81];% <=101 or 111

[Pts,Ims,nF] = fileload(dataset,inds);
X = Pts{1,1};
Y = Pts{1,2};
max_size = size(Ims{1},2);
order = randperm(30);
order = order(1:LX);

XX = X(order,:)/max_size;
YY = Y(1:LY,:)/max_size;

Pts{1,1} = X(order,:)';
Pts{1,2} = Y(1:LY,:)';

gt = zeros(LX,LY);
for i = 1:LX
    gt(i,order(i)) = 1;
end
asgT.alg = 'truth';
asgT.map = gt;asgT.X = gt;

%%
parKnl = st('alg', 'cmum_sc'); % type of affinity: only edge distance
[pars, algs] = gmPar(2);

parG = st('link', 'full'); % Delaunay triangulation for computing the graphs
gphs = newGphUs(Pts, parG);

[MM,~,~] = M_shape(Pts{1,1}',Pts{1,2}',1/8,2);
gphs{1}.sc = MM*1;

[KP, KQ] = conKnlGphPQU(gphs, parKnl);%%parKnl 可以用 pas 表示Pascal 数据集
K = conKnlGphKU(KP, KQ, gphs);
Ct = ones(size(KP));

asgTH.tim = 0;
asgTH.acc = 1;
asgTH.obj = gt(:)'*K*gt(:);

gphDs = gphU2Ds(gphs);
KQD = [KQ, KQ; KQ, KQ];


%% GA
asgGa = gm(K, Ct, asgT, pars{1}{:});

%% PM
asgPm = pm(K, KQ, gphs, asgT);

%% SM
asgSm = gm(K, Ct, asgT, pars{3}{:});

%% SMAC
asgSmac = gm(K, Ct, asgT, pars{4}{:});

%% IPFP-U
asgIpfpU = gm(K, Ct, asgT, pars{5}{:});

%% IPFP-S
asgIpfpS = gm(K, Ct, asgT, pars{6}{:});

%% RRWM
asgRrwm = gm(K, Ct, asgT, pars{7}{:});
%
%% FGM-U
asgFgmU = fgmU(KP, KQ, Ct, gphs, asgT, pars{8}{:});
%% FGM-D
asgFgmD = fgmD(KP, KQD, Ct, gphDs, asgT, pars{9}{:});
%%  FRGM
option.M_exist = 0.01;
option.alpha2 = 1;
option.lambda_EOT = 1/10;
option.niter_EOT = 20;
option.maxiter = 100;
option.active = 1;
option.q_norm = 1;
option.type = 'inner';

[M,~,~] = M_shape(XX,YY,1/8,2);
DX = M_points(XX,XX);
DY = M_points(YY,YY);
max_fact = max([max(DX(:)),max(DY(:))]);
DXX = DX/max_fact;DYY = DY/max_fact;
SX = 1./(DXX + 1*(DXX==0)).*(DXX > 0);
SX = SX/max(SX(:));

sx = neighbor1(XX);
sx = sx + sx';
sy = neighbor1(YY);
sy = sy + sy';
DXX_del = DXX.*sx;
DYY_del = DYY.*sy;

Map_ini = asgHun(-M);

sig = [0.5^2,0.5^2];
weight = [1,1];
[asgFRGM] = FRGM_Gen(Map_ini,M,DXX,DYY,SX,option,asgT,sig,weight);

%% MPM with/without delaunay
asgMPM = MPM_GM(X,Y,order,K,2500);


