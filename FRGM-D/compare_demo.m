%% comparison of ICP/CPD/JRMPC/cdfHC/GMM_reg/L2E/GLS/FGM-D

load('cpd_data2D_fish');

X = cpdfish{1,1};
LX = length(X(:,1));
[X,~] = normalize_point(X,1);

type = 'similar';% similar or nonrigid
switch type
    case 'similar'
        opt_rigid.noise = 1;
        opt_rigid.noise_sigma = 0.02;
        opt_rigid.noise_type = 'uniform';%'gaussian'
        opt_rigid.outlier = 100;
        opt_rigid.out_sigma = 0.2;
        opt_rigid.outlier_type = 'gaussian';%'uniform1''gaussian''gaussian1'
        
        theta = 0.1*pi;
        v = [cos(theta),sin(theta);
            -sin(theta),cos(theta)];
        t = 1*rand(1,2) + 0.5;
        s = 0.2;
        XT = rigid_affine_transform2D(X,v,t,s,opt_rigid);
        
        type_set = {'rigid2d','rigid','similar'};
        
    case 'nonrigid'
        sigma = 2;
        GX = nonrigid_ker(X,sigma,'rbf');
        theta = (0.0)*pi;
        R = [cos(theta),sin(theta);
            -sin(theta),cos(theta)];
        W = 0.5*randn(length(X(:,1)),2);
        t = 1*rand(1,2) + 0.5;
        scale = 0.5;
        
        opt_non.noise = 1;
        opt_non.noise_sigma = 0.02;
        opt_non.noise_type = 'uniform';%'gaussian'
        opt_non.outlier = max(floor(0.5*LX),1);
        opt_non.out_sigma = 0.25;
        opt_non.outlier_type = 'gaussian';%'uniform1''gaussian''gaussian1'
        
        XT = scale*nonrigid_kernel_trans(X,W,GX,R,opt_non);
        
        type_set = {'tps','nonrigid','nonrigid'};
end

figure,plot(XT(1:LX,1),XT(1:LX,2),'b.',X(:,1),X(:,2),'r.',XT((LX+1):end,1),XT((LX+1):end,2),'g.');

XT = normalize_point(XT,1);

if strcmp(type, 'similar')
    %% ICP
    XT = normalize_point(XT,1);
    X = normalize_point(X,1);
    [Ricp,Ticp,dataout]=icp2(XT', X',[],10,1);dataout = dataout';
    
    Xicp = (Ricp*X')'+repmat(Ticp',length(X(:,1)),1);
    
    plot_shapes(X,XT,dataout)

    [mseicp] = measurement(X,dataout,1:LX,[])
    
    %% JPRMPC
    XT = normalize_point(XT,1);
    X = normalize_point(X,1);
    [Rotation, TV] = JRMPC_2d(X,XT,100,1);
    
    [msejpr] = measurement([TV{1}(1:2,:)]',[TV{2}(1:2,:)]',1:LX,[])
end
%% cdfHC
XT = normalize_point(XT,1);
X = normalize_point(X,1);
cpdfish{1,1} = X;
cpdfish{1,2} = XT;
[T,En] = HC2Reg_TPS(cpdfish,1,200);

[msecdf] = measurement(T{1},T{2},1:LX,[])
%% L2E
XT = normalize_point(XT,1);
X = normalize_point(X,1);
order1 = 1:LX;
[idt,X_trans] = L2E_reg_ma(X(order1,:),XT,10,1,'iter');

[msel2e] = measurement(X_trans,XT,order1,[])
%% GLS
XT = normalize_point(XT,1);
X = normalize_point(X,1);
[Transform, C] = GLS_ma(X,XT,100,1,1);

[msegls] = measurement(Transform.Y,XT,1:LX,[])
%% GMMreg
[XT,Y_bound ]= normalize_point(XT,1);
[X,X_bound1 ]= normalize_point(X,1);
model = X;
scene = XT-repmat(mean(XT),length(XT(:,1)),1);
[config] = initialize_config(model,scene,type_set{1});
[param, transformed_model, history, config] = gmmreg_L2(config);

plot_shapes(model,scene,transformed_model)

[msegmm] = measurement(transformed_model,scene,1:LX,[])
%% CPD
opt.method=type_set{2};
opt.corresp=1;      % compute correspondence vector at the end of registration (not being estimated by default)
Transform=cpd_register(XT,X,opt);

[msegls] = measurement(Transform.Y,XT,1:LX,[])
%%
option.regist_trans = type_set{3};
option.regist_it = 50;
option.regist_rota = strcmp(type, 'similar');
option.regist_display = 1;
option.regist_normalize = 1;

option.GM_convex_or_non = [zeros(1,5),ones(1,option.regist_it)];%alternation of objective functions
option.GM_lambda1 = 1;% weight of unary potential
option.GM_lam_nonvex = 1;% weight of the non-convex pairwise potential
option.GM_lam_convex = 1;% weight of the convex pairwise potential
option.GM_initial = 'lap';% initialization 
option.GM_unary = 1;% linear unary potential
option.GM_connected = 'full';% fully connected graph
option.sigma = 4;% bandwidth of the RBF kernel
option.order = 1:LX;

[Map,para] = FRGM_pr(X,XT,option);

[msefrgmd] = measurement(para.X,XT,1:LX,[])

