% data load
data_name = 1;%1=beijing;2=whale;3=chinese;4=cpd_fish;5=fish_2d;
X = load_testdata(data_name);
LX = length(X(:,1));

% normalization of points
[X,X_bound] = normalize_point(X,1);

% generation of deformation, noise and outliers
sigma = 2;
GX = nonrigid_ker(X,sigma,'rbf');
theta = (0.0)*pi;
R = [cos(theta),sin(theta);
    -sin(theta),cos(theta)];
W = 0.5*randn(length(X(:,1)),2);
t = 1*rand(1,2) + 0.5;
%t=[0,0];
scale = 0.5;

opt_non.noise = 1;
opt_non.noise_sigma = 0.0;%0.02
opt_non.noise_type = 'gaussian';%'gaussian'
opt_non.outlier = 50;
opt_non.out_sigma = 0.25;
opt_non.outlier_type = 'gaussian';%'uniform1''gaussian''gaussian1'

P = [eye(LX,LX),zeros(LX,opt_non.outlier)];
Y = scale*nonrigid_kernel_trans(X,W,GX,R,opt_non);

figure,plot(Y(1:LX,1),Y(1:LX,2),'b.',X(:,1),X(:,2),'r.',Y((LX+1):end,1),Y((LX+1):end,2),'g.');

% randomly disrupt the order of points
LY = length(Y(:,1));
order = randperm(LX);
inlier_rate = 1;%in [0,1]: 1-inlier_rate of inliers are removed
order = order(1:min(floor(inlier_rate*LX)+1,LX));
%order = 1:LX;

X = X(order,:);

%% parameters

option.regist_trans = 'nonrigid';
option.regist_it = 50;%the max number of rigistration iterations
option.regist_rota = 0;% 1=: rotation invariant Shape Context; 0: otherwise
option.regist_display = 1;% 1:show the registrations; 0: otherwise
option.regist_normalize = 1;% 1: with normolization; 0: without
option.regist_save = 0;%save as gif or not

option.GM_convex_or_non = [zeros(1,5),ones(1,option.regist_it)];%alternation of objective functions
option.GM_lambda1 = 1;% weight of unary potential
option.GM_lam_nonvex = 10;% weight of the non-convex pairwise potential
option.GM_lam_convex = 1;% weight of the convex pairwise potential
option.GM_initial = 'lap';% initialization 
option.GM_unary = 1;% linear unary potential
option.GM_connected = 'full';% fully connected graph
option.order = order;

GXX = nonrigid_ker(X,4,'rbf');
tic;[Map,para] = FRGM_registration(X,Y,GXX,option);
toc;

figure,plot(1:length(para.regist_tol),para.regist_tol,'b.-','linewidth',2,'markersize',15);
hold on,plot(1:length(para.regist_err),para.regist_err,'r.-','linewidth',2,'markersize',15);
legend({'tol','err'},'FontSize',15)

