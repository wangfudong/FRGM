%% test cmu house/hotel
dataset = 1;% 1 for house (111 images), 2 for hotel (101 images)

LX = 20;%LX<=30;
LY = 30;%LY<=30;

inds = [1,51];% <=101 or 111

[Pts,Ims,nF] = fileload(dataset,inds);
X = Pts{1,1};
Y = Pts{1,2};
max_size = size(Ims{1},2);
order = randperm(30);
order = order(1:LX);

XX = X(order,:)/max_size;
YY = Y(1:LY,:)/max_size;

%% ATGM
option.Mexist = 0.1;     % weight of the unary term
option.maxiter_nonvex = 200; 
option.maxiter_convex = 100;
option.lam_nonvex = 10;  % weight of the nonconvex pairwise term
option.lam_sparse = 100; % weight of sparse_reuglarization term
option.lam_convex = 0.1; % weight of convex pairwise term 
option.geofunc = '1.21'; % the nonconvex objective function
option.unary = 1;        % the linear unary term
option.full = 1;         % 1 for fully-connected and 0 for delaunay triangulation

asgATGM = FRGM_E(XX,YY,order,option);
